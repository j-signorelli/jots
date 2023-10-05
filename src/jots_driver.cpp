#include "jots_driver.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

const string JOTSDriver::LINE = "-------------------------------------------------------------------";
const double JOTSDriver::TIME_TOLERANCE = 1e-14;

JOTSDriver::JOTSDriver(const Config& input, const int myid, const int num_procs, MPI_Comm in_comm)
: rank(myid),
  size(num_procs),
  comm(in_comm),
  sim(nullptr),
  adapter(nullptr),
  user_input(input),
  boundary_conditions(nullptr),
  initialized_bcs(false),
  ode_solver(nullptr),
  pmesh(nullptr),
  fe_coll(nullptr),
  fespace(nullptr),
  oper(nullptr),
  temp_T_gf(nullptr)
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
                                    
                                    

        )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.1" << endl << LINE << endl;
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
        ProcessTimeIntegration();
    //----------------------------------------------------------------------
    // Process preCICE (if using)
    if (user_input.UsingPrecice())
        ProcessPrecice();
    //----------------------------------------------------------------------
    // Process BoundaryConditions
    ProcessBoundaryConditions();
    //----------------------------------------------------------------------
    // Process LinearSolverSettings
    ProcessLinearSolverSettings();
    //----------------------------------------------------------------------
    // Process Output
    ProcessOutput();
    //---------------------------------------------------------------------
    // Create main solution vector from IC
    temp_T_gf->GetTrueDofs(T);
    //----------------------------------------------------------------------
    // Instantiate ConductionOperator, sending all necessary parameters
    if (rank == 0)
    {
        cout << LINE << endl;
        cout << "Initializing operator... ";
    }
    oper = new ConductionOperator(user_input, boundary_conditions, all_bdr_attr_markers, mat_props[MATERIAL_PROPERTY::DENSITY], mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT], mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY], *fespace, time);
    if (rank == 0)
        cout << "Done!" << endl;
    //----------------------------------------------------------------------
    // Instantiate OutputManager
    output = new OutputManager(rank, fespace, user_input, T, mat_props[MATERIAL_PROPERTY::DENSITY], mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT], mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);
}

void JOTSDriver::ProcessFiniteElementSetup()
{
    // Print appropriate solver type + get required material properties, add them to map
    if (rank == 0)
        cout << "Simulation Type: ";

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
        temp_T_gf = new ParGridFunction(fespace);
        
        //----------------------------------------------------------------------
        // Print initial condition info + set T_gf
        ConstantCoefficient coeff_IC(user_input.GetInitialTemp());
        temp_T_gf->ProjectCoefficient(coeff_IC);
        
        it_num = 0;
        time = 0.0;

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
        temp_T_gf = temp_visit_dc.GetParField("Temperature");

        // Take ownership of fe_coll and fespace (See MakeOwner docs)
        fe_coll = temp_T_gf->OwnFEC();
        fespace = temp_T_gf->ParFESpace();
        temp_T_gf->MakeOwner(NULL);

        // Now take ownership of mesh and temperature field
        // NOTE: If any further fields are added to restarts, they MUST be taken and deallocated here!!!!!
        // temp_T_gf deallocated below
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
}

void JOTSDriver::ProcessMaterialProperties()
{   
    if (rank == 0)
        cout << "\n";
    
    // Loop over INPUT material property info map from Config
    map<string, vector<string>> in_mat_props = user_input.GetMaterialPropertyInfoMap();
    map<string, vector<string>>::iterator it;

    for (it = in_mat_props.begin(); it != in_mat_props.end(); it++)
    {   
        // Print label
        string label = it->first;
        if (rank == 0)
            cout << label << ": ";
        
        // Get material property type
        MATERIAL_PROPERTY mp = Material_Property_Map.at(label);

        // Add material property to material property map
        vector<string> mp_info = it->second;
        vector<double> poly_coeffs;

        switch (Material_Model_Map.at(mp_info[0]))
        {
            case MATERIAL_MODEL::UNIFORM: // Uniform
                mat_props[mp] = new UniformProperty(stod(mp_info[1].c_str()));
                break;
            case MATERIAL_MODEL::POLYNOMIAL: // Polynomial

                for (int i = 1; i < mp_info.size(); i++)
                    poly_coeffs.push_back(stod(mp_info[i].c_str()));

                mat_props[mp] = new PolynomialProperty(poly_coeffs, *fespace);
                break;
            default:
                MFEM_ABORT("Unknown/Invalid material model specified");
                return;
        }

        // Print remaining portion of label
        cout << mat_props[mp]->GetInitString() << endl;
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
        cout << "Max Time: " << user_input.GetFinalTime() << endl;
    
    }
    //----------------------------------------------------------------------
    // Set ODE time integrator
    switch (Time_Scheme_Map.at(user_input.GetTimeSchemeLabel()))
    {
        case TIME_SCHEME::EULER_IMPLICIT:
            ode_solver = new BackwardEulerSolver;
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            ode_solver = new ForwardEulerSolver;
            break;
        case TIME_SCHEME::RK4:
            ode_solver = new RK4Solver;
    }
    //----------------------------------------------------------------------
    // Set time stuff
    dt = user_input.Getdt();
    tf = user_input.GetFinalTime();
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

    // Confirm user input matches mesh bdr_attributes...
    // Check count
    if (pmesh->bdr_attributes.Size() != user_input.GetBCCount())
    {
        if (rank == 0)
            MFEM_ABORT("Input file BC count and mesh BC counts do not match.");
        return;
    }
    
    // Check one-to-oneness
    // Also create an array of indices that match mesh bdr attributes to input bdr attributes
    // ^ in event that user inputs bdr attributes in different order than mesh
    vector<int> bdr_index;
    for (size_t i = 0; i < user_input.GetBCCount(); i++)
    {
        int attr = pmesh->bdr_attributes[i];
        bool one_to_one = false;

        for (size_t j = 0; j < user_input.GetBCCount(); j++)
        {   
            pair<int, vector<string>> bc = user_input.GetBCInfo(j);
            if (attr == bc.first)
            {
                one_to_one = true;
                bdr_index.push_back(j);
                j = user_input.GetBCCount();
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
        pair<int, vector<string>> bc = user_input.GetBCInfo(bdr_index[i]);
        double value;
        string mesh_name;
        double amp;
        double afreq;
        double phase;
        double shift;
        switch (Boundary_Condition_Map.at(bc.second[0]))
        {
            case BOUNDARY_CONDITION::ISOTHERMAL:
                value = stod(bc.second[1].c_str());
                boundary_conditions[i] =  new UniformConstantIsothermalBC(bc.first, value);
                break;
            case BOUNDARY_CONDITION::HEATFLUX:
                value = stod(bc.second[1].c_str());
                boundary_conditions[i] =  new UniformConstantHeatFluxBC(bc.first, value);
                break;
            case BOUNDARY_CONDITION::PRECICE_ISOTHERMAL:
                mesh_name = bc.second[1];
                value = stod(bc.second[2].c_str());
                boundary_conditions[i] =  new PreciceIsothermalBC(bc.first, *fespace, mesh_name, value);
                precice_bc_indices.push_back(i);
                break;
            case BOUNDARY_CONDITION::PRECICE_HEATFLUX:
                mesh_name = bc.second[1];
                value = stod(bc.second[2].c_str());
                boundary_conditions[i] =  new PreciceHeatFluxBC(bc.first, *fespace, mesh_name, value);
                precice_bc_indices.push_back(i);
                break;
            case BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL:
                amp = stod(bc.second[1].c_str());
                afreq = stod(bc.second[2].c_str());
                phase = stod(bc.second[3].c_str());
                shift = stod(bc.second[4].c_str());
                boundary_conditions[i] = new UniformSinusoidalIsothermalBC(bc.first, time, amp, afreq, phase, shift);
                break;
            case BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX:
                amp = stod(bc.second[1].c_str());
                afreq = stod(bc.second[2].c_str());
                phase = stod(bc.second[3].c_str());
                shift = stod(bc.second[4].c_str());
                boundary_conditions[i] = new UniformSinusoidalHeatFluxBC(bc.first, time, amp, afreq, phase, shift);
                break;
            default:
                MFEM_ABORT("Invalid/Unknown boundary condition specified: '" + bc.second[0] + "'");
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
    for (size_t i = 0; i < user_input.GetBCCount(); i++)
    {
        Array<int> bdr_attr(user_input.GetBCCount());
        bdr_attr = 0;
        bdr_attr[i] = 1;
        all_bdr_attr_markers[i] = bdr_attr;
   }
    //----------------------------------------------------------------------
    // Print BCs
    for (size_t i = 0; i < user_input.GetBCCount(); i++)
    {   
        BoundaryCondition* bc = boundary_conditions[i];

        if (rank == 0)
        {
            cout << "Boundary Attribute " << bc->GetBdrAttr() << ": " << bc->GetInitString() << endl;
        }
    }
}

void JOTSDriver::ProcessLinearSolverSettings()
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

void JOTSDriver::ProcessOutput()
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
    // Initialize material properties (only pertinent for non-uniform non-constant ones)
    if (rank == 0)
        cout << "Initializing material properties field...";

    UpdateMatProps();

    if (rank == 0)
        cout << " Done!" << endl;
    

    // Initialize the ODE Solver
    if (rank == 0)
        cout << "Initializing solver...";

    ode_solver->Init(*oper);

    if (rank == 0)
        cout << " Done!" << endl;



    double precice_dt = 0;
    // Implicit coupling:
    double precice_saved_time = 0;
    double precice_saved_it_num = 0;

    // Output IC as ParaView
    if (it_num == 0)
        output->WriteVizOutput(it_num, time);


    // If using precice: initialize + get/send initial data if needed
    // Make actual precice interface calls directly here for readability/clarity
    if (user_input.UsingPrecice())
    {
        precice_dt = adapter->Interface()->initialize();
        if (adapter->Interface()->isActionRequired(PreciceAdapter::cowid))
        {
            adapter->WriteData(T, mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);
            adapter->Interface()->markActionFulfilled(PreciceAdapter::cowid);
        }
        adapter->Interface()->initializeData();
    }
    
    bool viz_out = false;
    bool res_out = false;

    while ( (!user_input.UsingPrecice() && time < tf) 
        || (user_input.UsingPrecice() && adapter->Interface()->isCouplingOngoing()))//Main Solver Loop - use short-circuiting
    {

        if (user_input.UsingPrecice())
        {
            if (adapter->Interface()->isActionRequired(PreciceAdapter::cowic))
            {
                precice_saved_time = time;
                precice_saved_it_num = it_num;
                adapter->SaveOldState(T);
                adapter->Interface()->markActionFulfilled(PreciceAdapter::cowic);
            }
            if (adapter->Interface()->isReadDataAvailable())
                adapter->GetReadData();
        }


        // Update BCs and update BLFs + LFs in operator, if required
        PreprocessIteration();


        // Update timestep if needed
        if (user_input.UsingPrecice())
        {
            if (precice_dt < dt)
                dt = precice_dt;
        }

        // Step in time - time automatically updated
        // NOTE: Do NOT use ANY ODE-Solvers that update dt
        ode_solver->Step(T, time, dt);        
        it_num++; // increment it_num

        // Leave solver if time now greater than tf -- solution is done
        if (time > tf && abs(time-tf) > TIME_TOLERANCE)
            continue;
        
        // Update material properties with new temperature
        UpdateMatProps();

        // Write any preCICE data, reload state if needed
        if (user_input.UsingPrecice())
        {
            if (adapter->Interface()->isWriteDataRequired(dt))
                adapter->WriteData(T, mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);

            // Advance preCICE
            adapter->Interface()->advance(dt);


            // Implicit coupling
            if (adapter->Interface()->isActionRequired(PreciceAdapter::coric))
            {
                time = precice_saved_time;
                it_num = precice_saved_it_num;
                adapter->ReloadOldState(T);
                UpdateMatProps(); // Update material props w/ reloaded field
                adapter->Interface()->markActionFulfilled(PreciceAdapter::coric);
                continue; // skip printing of timestep info AND outputting
            }
        }
   

        // Print current timestep information:
        if (rank == 0)
        {
            printf("Step #%10i || Time: %1.6e out of %-1.6e || dt: %1.6e \n", it_num, time, tf, dt);
                
        }
        
        // Check if blow up
        if (T.Max() > 1e10)
        {
            MFEM_ABORT("JOTS has blown up");
            return;
        }

        // Output
        viz_out = user_input.GetVisFreq() != 0 && it_num % user_input.GetVisFreq() == 0;
        res_out = user_input.GetRestartFreq() != 0 && it_num % user_input.GetRestartFreq() == 0;
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
    
    // Finalize preCICE if used
    if (user_input.UsingPrecice())
        adapter->Interface()->finalize();


}

void JOTSDriver::UpdateMatProps()
{
    if (!mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]->IsConstant())
        mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]->UpdateCoeff(T); // Update k coefficient

    if (!mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT]->IsConstant())
        mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT]->UpdateCoeff(T); // Update C coefficient
}

void JOTSDriver::PreprocessIteration()
{
    // If specific heat not constant, it may have changed -> Update mass BLF
    if (!mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT]->IsConstant())
        oper->UpdateMass();

    // If conductivity not constant, it may have changed -> Update stiffness BLF
    if (!mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]->IsConstant())
        oper->UpdateStiffness();

    // Update time-dependent/non-constant boundary condition coefficients
    temp_T_gf->SetFromTrueDofs(T);

    bool n_changed = false;
    bool d_changed = false;
    for (size_t i = 0; i < user_input.GetBCCount(); i++)
    {   
        // Update coefficients (could be preCICE calls, could be SetTime calls, etc.)
        if (!boundary_conditions[i]->IsConstant() || !initialized_bcs) // If not constant in time or not yet initialized for first iteration
        {
            boundary_conditions[i]->UpdateCoeff();

            if (boundary_conditions[i]->IsEssential())
            {
                // Project correct values on boundary for essential BCs
                d_changed = true;
                temp_T_gf->ProjectBdrCoefficient(boundary_conditions[i]->GetCoeffRef(), all_bdr_attr_markers[i]);
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
        // Apply changed Dirichlet BCs to T
        temp_T_gf->GetTrueDofs(T);
        // Here is where can set tmp_du_dt to be used if HO wanted
        // For higher-order, need T from n-1 timestep in addition to n timestep? is this worth doing?
        // Need to save previous timestep temperature in restarts
    }

    // If any non-constant Neumann terms, update Neumann linear form
    if (n_changed)
        oper->UpdateNeumannTerm();

}

JOTSDriver::~JOTSDriver()
{
    delete sim;
    delete adapter;
    for (map<MATERIAL_PROPERTY, MaterialProperty*>::iterator it = mat_props.begin(); it != mat_props.end(); it++)
        delete it->second;
    for (size_t i = 0; i < user_input.GetBCCount(); i++)
        delete boundary_conditions[i];
    delete[] boundary_conditions;
    delete[] all_bdr_attr_markers;
    delete ode_solver;
    delete pmesh;
    delete fe_coll;
    delete fespace;
    delete oper;
    delete output;
    delete temp_T_gf;

}