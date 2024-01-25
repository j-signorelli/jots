#include "jots_driver.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

const string JOTSDriver::LINE = "-------------------------------------------------------------------";

JOTSDriver::JOTSDriver(const Config& input, const int myid, const int num_procs, MPI_Comm in_comm)
: rank(myid),
  size(num_procs),
  comm(in_comm),
  it_num(0),
  time(0.0),
  dt(0.0),
  max_timesteps(0),
  jots_iterator(nullptr),
  adapter(nullptr),
  user_input(input),
  boundary_conditions(nullptr),
  all_bdr_attr_markers(nullptr),
  initialized_bcs(false),
  mat_props(nullptr),
  pmesh(nullptr),
  fe_coll(nullptr),
  fespace(nullptr),
  u_0_gf(nullptr)
{   
    if (rank == 0)
    {
        cout << LINE << endl;
        cout << R"(

              _  ____ _______ _____ 
             | |/ __ \__   __/ ____|
             | | |  | | | | | (___  
         _   | | |  | | | |  \___ \
        | |__| | |__| | | |  ____) |
         \____/ \____/  |_| |_____/ 
                                    
                                    

        )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.2" << endl << LINE << endl;
    }
    //----------------------------------------------------------------------
    // Print config file
    if (rank == 0)
      cout << "Configuration file: " << user_input.GetInputFile() << endl;
    //----------------------------------------------------------------------
    // Process FiniteElementSetup
    ProcessFiniteElementSetup();
    //----------------------------------------------------------------------
    // Process MaterialProperties in config file
    ProcessMaterialProperties();
    //----------------------------------------------------------------------
    // Process TimeIntegration (if using (required for unsteady))
    if (user_input.UsingTimeIntegration())
    {
        ProcessTimeIntegration();
    }
    else // otherwise set cycle to -2, max timesteps to -1
    {
        it_num = -2;
        max_timesteps = -1;
        // Ensures output appropriate for time-independent runs, as MFEM requires cycle = -1
        // Time-independent runs only complete inner iterations w/ a single outer iteration
        // As opposed to time-dependent runs having outer-iterations being a timestep
    }
    //----------------------------------------------------------------------
    // Process preCICE (if using)
    if (user_input.UsingPrecice())
        ProcessPrecice();
    //----------------------------------------------------------------------
    // Process BoundaryConditions
    ProcessBoundaryConditions();
    //----------------------------------------------------------------------
    // Print LinearSolverSettings
    PrintLinearSolverSettings();
    //----------------------------------------------------------------------
    // Print NewtonSolverSettings (if using)
    if (user_input.UsingNewton())
        PrintNewtonSolverSettings();
    
    //----------------------------------------------------------------------
    // Print Output settings (if unsteady)
    if (user_input.UsingTimeIntegration())
        PrintOutput();
    //---------------------------------------------------------------------
    // Create JOTSIterator object/s
    if (rank == 0)
    {
        cout << LINE << endl;
        cout << "Initializing iterator... ";
    }
    switch (Simulation_Type_Map.at(user_input.GetSimTypeLabel()))
    {
        case SIMULATION_TYPE::UNSTEADY:
            jots_iterator = new ConductionOperator(user_input,
                                             boundary_conditions,
                                             all_bdr_attr_markers,
                                             *mat_props[MATERIAL_PROPERTY::DENSITY],
                                             *mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT],
                                             *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY],
                                             *fespace,
                                             time,
                                             dt);
            break;
        case SIMULATION_TYPE::STEADY:
            jots_iterator = new SteadyConductionOperator(user_input, 
                                                    boundary_conditions,
                                                    all_bdr_attr_markers,
                                                    *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY],
                                                    *fespace);
            break;
    }
    //----------------------------------------------------------------------
    if (rank == 0)
        cout << "Done!" << endl;
    //----------------------------------------------------------------------

    // Initialize material properties w/ solution field (only pertinent for non-uniform non-constant ones)
    if (rank == 0)
        cout << "Initializing material properties field...";

    UpdateAndApplyMatProps();

    if (rank == 0)
        cout << " Done!" << endl;
}

void JOTSDriver::ProcessFiniteElementSetup()
{
    // Print appropriate solver type + get required material properties, add them to map
    if (rank == 0)
        cout << "Simulation Type: " << user_input.GetSimTypeLabel() << endl;

    // If not restart, refine mesh and initialize; else load VisItDataCollection
    if (!user_input.UsesRestart())
    {
        if (rank == 0)
            cout << "Non-restart simulation..." << endl;
        //----------------------------------------------------------------------
        // Create serial mesh
        Mesh* mesh = new Mesh(user_input.GetMeshFile().c_str(), 1);//pass one to generate edges
        dim = mesh->Dimension();
        //----------------------------------------------------------------------
        // Print mesh info
        if (rank == 0)
        {
            cout << "\n";
            cout << "Mesh File: " << user_input.GetMeshFile() << endl;
        }
        //----------------------------------------------------------------------
        // Refine mesh in serial
        int ser_ref = user_input.GetSerialRefine();
        for (int lev = 0; lev < ser_ref; lev++)
        {
            mesh->UniformRefinement();
        }
        //----------------------------------------------------------------------
        // Print serial refinements
        if (rank == 0)
            cout << "Serial refinements of mesh completed: " << ser_ref << endl;
        //----------------------------------------------------------------------
        // Complete parallel decomposition of serial mesh + refine parallel
        pmesh = new ParMesh(comm, *mesh);
        delete mesh;
        int par_ref = user_input.GetParallelRefine();

        for (int lev = 0; lev < par_ref; lev++)
        {
            pmesh->UniformRefinement();
        }
        //----------------------------------------------------------------------
        // Print parallel refinements
        if (rank == 0)
        {
            cout << "Parallel refinements of mesh completed: " << par_ref << endl;
        }
        //----------------------------------------------------------------------
        // Define parallel FE space on parallel mesh
        fe_coll = new H1_FECollection(user_input.GetFEOrder(), dim);
        
        fespace = new ParFiniteElementSpace(pmesh, fe_coll);
        
        //----------------------------------------------------------------------
        // Create solution GF
        u_0_gf = new ParGridFunction(fespace);
        
        //----------------------------------------------------------------------
        // Print initial condition info + set T_gf
        ConstantCoefficient coeff_IC(user_input.GetInitialTemp());
        u_0_gf->ProjectCoefficient(coeff_IC);
        

        if (rank == 0)
            cout << "Initial temperature field: " << user_input.GetInitialTemp() << endl;

    }
    //----------------------------------------------------------------------
    // Else using restart:
    else
    {   
        if (rank == 0)
            cout << "Restart simulation..." << endl;
        //----------------------------------------------------------------------
        // Read in VisItDataCollection
        VisItDataCollection temp_visit_dc(comm, user_input.GetRestartPrefix());
        temp_visit_dc.Load(user_input.GetRestartCycle());

        if (temp_visit_dc.Error())
            MFEM_ABORT("Unable to load restart data collection");
        
        //----------------------------------------------------------------------
        // Get the parallel mesh
        pmesh = dynamic_cast<ParMesh*>(temp_visit_dc.GetMesh());
        dim = pmesh->Dimension();

        //----------------------------------------------------------------------
        // Get the temperature grid function
        u_0_gf = temp_visit_dc.GetParField("Temperature");

        // Take ownership of fe_coll and fespace (See MakeOwner docs)
        fe_coll = u_0_gf->OwnFEC();
        fespace = u_0_gf->ParFESpace();
        u_0_gf->MakeOwner(NULL);

        // Now take ownership of mesh and temperature field
        // NOTE: If any further fields are added to restarts, they MUST be taken and deallocated here!!!!!
        // pmesh deallocated in destructor
        
        temp_visit_dc.SetOwnData(false);

        // Set the it_num and time
        it_num = temp_visit_dc.GetCycle();
        time = temp_visit_dc.GetTime();


        if (rank == 0)
        {
            cout << "Input Restart Root: " << user_input.GetRestartPrefix() << "_" << user_input.GetRestartCycle() << endl;
            cout << "\tTime: " << time << endl;
            cout << "\tCycle: " << it_num << endl;
        }
    }
    //----------------------------------------------------------------------
    // Print mesh attributes
    if (rank == 0)
    {
        cout << "Problem Dimension: " << dim << endl;
        cout << "Mesh Boundary Attributes: <";
        for (int i = 0; i < pmesh->bdr_attributes.Size(); i++)
        {
            if (i != 0)
                cout << ",";
            cout << pmesh->bdr_attributes[i];
        }

        cout << ">" << endl;
    }
    //----------------------------------------------------------------------
    // Print number of unknowns
    HYPRE_BigInt fe_size = fespace->GlobalTrueVSize();
    if (rank == 0)
        cout << "Number of temperature nodes: " << fe_size << endl;

    //----------------------------------------------------------------------
    // Instantiate OutputManager
    output = new OutputManager(rank, *fespace, user_input);

    //----------------------------------------------------------------------
    // Set solution vector from IC
    u_0_gf->GetTrueDofs(u);

    // Add solution vector to OutputManager
    output->RegisterSolutionVector("Temperature", u);

}

void JOTSDriver::ProcessMaterialProperties()
{   
    if (rank == 0)
        cout << "\n";

    // Instantiate material property array w/ max number of mat props
    mat_props = new MaterialProperty*[Material_Property_Map.size()];
    for (size_t i = 0; i < Material_Property_Map.size(); i++)
        mat_props[i] = nullptr;

    // Loop over INPUT material property info map from Config
    // Get keys
    vector<string> in_mat_prop_labels = user_input.GetMaterialPropertyKeys();

    for (size_t i = 0; i < in_mat_prop_labels.size(); i++)
    {   
        // Print label
        string label = in_mat_prop_labels[i];
        if (rank == 0)
            cout << label << ": ";
        
        // Get material property type
        MATERIAL_PROPERTY mp = Material_Property_Map.at(label);

        // Add material property to material property array
        vector<string> mp_info = user_input.GetMaterialPropertyInfo(label);
        vector<double> poly_coeffs;

        MATERIAL_MODEL mm = Material_Model_Map.at(mp_info[0]);
        switch (mm)
        {
            case MATERIAL_MODEL::UNIFORM: // Uniform
                mat_props[mp] = new UniformProperty(stod(mp_info[1].c_str()));
                break;
            case MATERIAL_MODEL::POLYNOMIAL: // Polynomial

                for (size_t i = 1; i < mp_info.size(); i++)
                    poly_coeffs.push_back(stod(mp_info[i].c_str()));

                mat_props[mp] = new PolynomialProperty(poly_coeffs, *fespace);
                break;
            default:
                MFEM_ABORT("Unknown/Invalid material model specified");
                return;
        }

        // Print remaining portion of label
        if (rank == 0)
            cout << mat_props[mp]->GetInitString() << endl;

        // Register material property output
        output->RegisterCoefficient(label, mat_props[mp]->GetCoeffRef());
    }

}

void JOTSDriver::ProcessTimeIntegration()
{   
    //----------------------------------------------------------------------
    // Print time integration information
    if (rank == 0)
    {
        cout << "\n";
        cout << "Time Scheme: " << user_input.GetTimeSchemeLabel() << endl;
        cout << "Time Step: " << user_input.Getdt() << endl;
        cout << "Max Timesteps: " << user_input.GetMaxTimesteps() << endl;
    
    }
    //----------------------------------------------------------------------
    // Set dt and max timesteps
    dt = user_input.Getdt();
    max_timesteps = user_input.GetMaxTimesteps();
}

void JOTSDriver::ProcessPrecice()
{
    // Print any precice info + instantiate adapter object if needed
    if (rank == 0)
    {
        cout << "\n";
        cout << "Using preCICE!" << endl;
        cout << "preCICE Participant Name: " << user_input.GetPreciceParticipantName() << endl;
        cout << "preCICE Config File: " << user_input.GetPreciceConfigFile() << endl;
    }

    adapter = new PreciceAdapter(user_input.GetPreciceParticipantName(), user_input.GetPreciceConfigFile(), rank, size, comm);

    if (adapter->GetDimension() != dim)
    {
        MFEM_ABORT("preCICE dimensions and mesh file dimensions are not the same!");
        return;
    }
}

void JOTSDriver::ProcessBoundaryConditions()
{
    // Verify input config BCs appropriately match input mesh BCs
    if (rank == 0)
        cout << "\n";

    // Confirm user input matches mesh bdr_attributes size...
    // Check count
    if (pmesh->bdr_attributes.Size() != user_input.GetBCCount())
    {
        if (rank == 0)
            MFEM_ABORT("Input file BC count and mesh BC counts do not match.");
        return;
    }
    
    // Check one-to-oneness
    vector<int> bc_keys = user_input.GetBCKeys();
    for (int i = 0; i < pmesh->bdr_attributes.Size(); i++)
    {
        int attr = pmesh->bdr_attributes[i];
        bool one_to_one = false;

        for (size_t j = 0; j < bc_keys.size(); j++)
        {
            if (attr == bc_keys[j])
            {
                one_to_one = true;
                j = bc_keys.size();
            }
        }
        if (!one_to_one)
        {   
            if (rank == 0)
            {
                stringstream sstm;
                sstm << "No matching boundary attribute in config file for attribute " << attr;
                MFEM_ABORT(sstm.str());
            }
            return;
        }
    }
    //----------------------------------------------------------------------
    // Instantiate boundary conditions + save any precice bc indices
    vector<int> precice_bc_indices;
    boundary_conditions = new BoundaryCondition*[user_input.GetBCCount()];
    for (int i = 0; i < user_input.GetBCCount(); i++)
    {   
        int attr = pmesh->bdr_attributes[i];
        vector<string> bc_info = user_input.GetBCInfo(attr);
        double value;
        string mesh_name;
        double amp;
        double afreq;
        double phase;
        double shift;
        switch (Boundary_Condition_Map.at(bc_info[0]))
        {
            case BOUNDARY_CONDITION::ISOTHERMAL:
                value = stod(bc_info[1].c_str());
                boundary_conditions[i] =  new UniformConstantIsothermalBC(attr, value);
                break;
            case BOUNDARY_CONDITION::HEATFLUX:
                value = stod(bc_info[1].c_str());
                boundary_conditions[i] =  new UniformConstantHeatFluxBC(attr, value);
                break;
            case BOUNDARY_CONDITION::PRECICE_ISOTHERMAL:
                mesh_name = bc_info[1];
                value = stod(bc_info[2].c_str());
                boundary_conditions[i] =  new PreciceIsothermalBC(attr, *fespace, mesh_name, value);
                precice_bc_indices.push_back(i);
                break;
            case BOUNDARY_CONDITION::PRECICE_HEATFLUX:
                mesh_name = bc_info[1];
                value = stod(bc_info[2].c_str());
                boundary_conditions[i] =  new PreciceHeatFluxBC(attr, *fespace, mesh_name, value);
                precice_bc_indices.push_back(i);
                break;
            case BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL:
                amp = stod(bc_info[1].c_str());
                afreq = stod(bc_info[2].c_str());
                phase = stod(bc_info[3].c_str());
                shift = stod(bc_info[4].c_str());
                boundary_conditions[i] = new UniformSinusoidalIsothermalBC(attr, amp, afreq, phase, shift);
                break;
            case BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX:
                amp = stod(bc_info[1].c_str());
                afreq = stod(bc_info[2].c_str());
                phase = stod(bc_info[3].c_str());
                shift = stod(bc_info[4].c_str());
                boundary_conditions[i] = new UniformSinusoidalHeatFluxBC(attr, amp, afreq, phase, shift);
                break;
            default:
                MFEM_ABORT("Invalid/Unknown boundary condition specified: '" + bc_info[0] + "'");
                return;
        }
    }
    //----------------------------------------------------------------------
    // Send precice bcs to adapter
    if (user_input.UsingPrecice())
        adapter->AddPreciceBCs(boundary_conditions, precice_bc_indices);
    //----------------------------------------------------------------------
    // Prepare BC attr arrays for applying coefficients
    all_bdr_attr_markers = new Array<int>[user_input.GetBCCount()];
    for (int i = 0; i < user_input.GetBCCount(); i++)
    {
        Array<int> bdr_attr(user_input.GetBCCount());
        bdr_attr = 0;
        bdr_attr[i] = 1;
        all_bdr_attr_markers[i] = bdr_attr;
    }
    //----------------------------------------------------------------------
    // Print BCs
    for (int i = 0; i < user_input.GetBCCount(); i++)
    {   
        BoundaryCondition* bc = boundary_conditions[i];

        if (rank == 0)
        {
            cout << "Boundary Attribute " << bc->GetBdrAttr() << ": " << bc->GetInitString() << endl;
        }
    }
}

void JOTSDriver::PrintLinearSolverSettings()
{
    // Print linear solver settings
    if (rank == 0)
    {
        cout << "\n";
        cout << "Linear Solver: " << user_input.GetSolverLabel() << endl;
        cout << "Preconditioner: " << user_input.GetPrecLabel() << endl;
        cout << "Max Iterations: " << user_input.GetMaxIter() << endl;
        cout << "Absolute Tolerance: " << user_input.GetAbsTol() << endl;
        cout << "Relative Tolerance: " << user_input.GetRelTol() << endl;
    }
}

void JOTSDriver::PrintNewtonSolverSettings()
{
    // Print Newton solver settings
    if (rank == 0)
    {
        cout << "\n";
        cout << "Newton Solver" << endl;
        cout << "Max Iterations: " << user_input.GetNewtonMaxIter() << endl;
        cout << "Absolute Tolerance: " << user_input.GetNewtonAbsTol() << endl;
        cout << "Relative Tolerance: " << user_input.GetNewtonRelTol() << endl;
    }
}

void JOTSDriver::PrintOutput()
{
    // Print output settings
    if (rank == 0)
    {
        cout << "\n\n";
        cout << "Restart Frequency: " << user_input.GetRestartFreq() << endl;
        cout << "Visualization Frequency: " << user_input.GetVisFreq() << endl;
    }
}

void JOTSDriver::Run()
{   

    // Output IC to ParaView if time-integrating simulation
    if (it_num == 0 && user_input.UsingTimeIntegration())
        output->WriteVizOutput(it_num, time);

    // If using precice: initialize + get/send initial data if needed
    // Make actual precice interface calls directly here for readability/clarity
    if (user_input.UsingPrecice())
    {
        precice_dt = adapter->Interface()->initialize();
        if (adapter->Interface()->isActionRequired(PreciceAdapter::cowid))
        {
            adapter->WriteData(u, mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);
            adapter->Interface()->markActionFulfilled(PreciceAdapter::cowid);
        }
        adapter->Interface()->initializeData();
    }
    
    while ((!user_input.UsingPrecice() && it_num < max_timesteps) 
    || (user_input.UsingPrecice() && adapter->Interface()->isCouplingOngoing())) // use short-circuiting
    {
        // Handle preCICE calls
        if (user_input.UsingPrecice())
        {
            // Implicit coupling: save state
            if (adapter->Interface()->isActionRequired(PreciceAdapter::cowic))
            {
                precice_saved_time = time;
                precice_saved_it_num = it_num;
                adapter->SaveOldState(u);
                adapter->Interface()->markActionFulfilled(PreciceAdapter::cowic);
            }

            // Get read data
            if (adapter->Interface()->isReadDataAvailable())
                adapter->GetReadData();

            // Update timestep if needed
            if (precice_dt < dt)
                dt = precice_dt;
        }

        // Update and Apply BCs
        UpdateAndApplyBCs();

        // Iterate
        Iteration();

        // Write any preCICE data, reload state if needed
        if (user_input.UsingPrecice())
        {
            if (adapter->Interface()->isWriteDataRequired(dt))
                adapter->WriteData(u, mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);

            // Advance preCICE
            adapter->Interface()->advance(dt);


            // Implicit coupling
            if (adapter->Interface()->isActionRequired(PreciceAdapter::coric))
            {
                time = precice_saved_time;
                it_num = precice_saved_it_num;
                adapter->ReloadOldState(u);
                UpdateAndApplyMatProps(); // Update material props w/ reloaded field
                adapter->Interface()->markActionFulfilled(PreciceAdapter::coric);
                continue; // skip printing of timestep info AND outputting
            }
        }
   
        PostprocessIteration();
    }
    
    // Finalize preCICE if used
    if (user_input.UsingPrecice())
        adapter->Interface()->finalize();


}

void JOTSDriver::UpdateAndApplyMatProps()
{
    for (size_t i = 0; i < Material_Property_Map.size(); i++)
    {
        MATERIAL_PROPERTY mp = MATERIAL_PROPERTY(i);
        if (mat_props[mp] != nullptr && !mat_props[mp]->IsConstant())
        {
            // Update coefficient using current solution field, ie: k=k(T)
            mat_props[mp]->UpdateCoeff(u);
            
            // Update any BLFs affected by changed coefficient (Apply)
            jots_iterator->ProcessMatPropUpdate(mp);
        }
    }
}

void JOTSDriver::UpdateAndApplyBCs()
{
    // Update time-dependent/non-constant boundary condition coefficients
    // Get GF from current solution
    u_0_gf->SetFromTrueDofs(u);

    bool n_changed = false;
    bool d_changed = false;
    for (int i = 0; i < user_input.GetBCCount(); i++)
    {   
        // Update coefficients (could be preCICE calls, could be SetTime calls, etc.)
        if (!boundary_conditions[i]->IsConstant() || !initialized_bcs) // If not constant in time or not yet initialized for first iteration
        {
            boundary_conditions[i]->UpdateCoeff(time);

            if (boundary_conditions[i]->IsEssential())
            {
                // Project correct values on boundary for essential BCs
                d_changed = true;
                u_0_gf->ProjectBdrCoefficient(boundary_conditions[i]->GetCoeffRef(), all_bdr_attr_markers[i]);
            }
            else
            {
                n_changed = true; // Flag that Neumann LF must be updated
            }
        }
    }

    if (!initialized_bcs)
        initialized_bcs = true;

    if (d_changed)
    {
        // Apply changed Dirichlet BCs to u
        u_0_gf->GetTrueDofs(u);
        // Here is where can set tmp_du_dt to be used if HO wanted
        // For higher-order, need T from n-1 timestep in addition to n timestep? is this worth doing?
        // Need to save previous timestep temperature in restarts
    }

    // If any non-constant Neumann terms, update Neumann linear form
    if (n_changed)
        jots_iterator->UpdateNeumann();
}


void JOTSDriver::Iteration()
{
    // Step simulation
    // NOTE: Do NOT use ANY ODE-Solvers that update dt
    jots_iterator->Iterate(u);
    it_num++; // increment it_num - universal metric
    
    // Update material properties with new solution field - important to do before preCICE writes
    UpdateAndApplyMatProps();
}

void JOTSDriver::PostprocessIteration()
{   

    // Print current timestep information if using TimeIntegration:
    if (rank == 0 && user_input.UsingTimeIntegration())
    {
        printf("Step #%10i || Time: %1.6e out of %-1.6e || dt: %1.6e \n", it_num, time, dt*max_timesteps, dt);
    }
    
    // Check if blow up
    if (u.Max() > 1e10)
    {
        MFEM_ABORT("JOTS has blown up");
        return;
    }

    // Write any output files if required
    // If unsteady simulation, then output at input frequencies
    // If steady, then output ONLY IF viz_freq and restart_freq != 0
    if (!user_input.UsingTimeIntegration())
    {
        if (user_input.GetVisFreq() != 0)
        {
            if (rank == 0)
            {
                cout << LINE << endl;
                cout << "Saving Paraview Data..." << endl;
            }
            output->WriteVizOutput(it_num, time);
        }
        if (user_input.GetRestartFreq() != 0)
        {
            if (rank == 0)
            {
                cout << "Saving Restart File..." << endl;
                cout << LINE << endl;
            }
            output->WriteRestartOutput(it_num, time);
        }
    }
    else
    {
        bool viz_out = user_input.GetVisFreq() != 0 && it_num % user_input.GetVisFreq() == 0;
        bool res_out = user_input.GetRestartFreq() != 0 && it_num % user_input.GetRestartFreq() == 0;
        if (viz_out || res_out)
        {
            if (rank == 0)
                cout << LINE << endl;
                
            if (viz_out)
            {
                if (rank == 0)
                    cout << "Saving Paraview Data: Cycle " << it_num << endl;
                output->WriteVizOutput(it_num, time);
            }

            if (res_out)
            {
                if (rank == 0)
                    cout << "Saving Restart File: Cycle " << it_num << endl;
                output->WriteRestartOutput(it_num, time);
            }

            if (rank == 0)
                cout << LINE << endl;
        }
    }
}

JOTSDriver::~JOTSDriver()
{
    delete jots_iterator;
    delete adapter;
    for (size_t i = 0; i < Material_Property_Map.size(); i++)
        delete mat_props[MATERIAL_PROPERTY(i)];
    delete[] mat_props;
    for (int i = 0; i < user_input.GetBCCount(); i++)
        delete boundary_conditions[i];
    delete[] boundary_conditions;
    delete[] all_bdr_attr_markers;
    delete pmesh;
    delete fe_coll;
    delete fespace;
    delete output;
    delete u_0_gf;

}