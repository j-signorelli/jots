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
  jots_iterator(), // initialize all ptrs to nullptr
  u(),
  precice_interface(nullptr),
  user_input(input),
  boundary_conditions(),
  all_bdr_attr_markers(nullptr),
  initialized_bcs(false),
  mat_props(),
  pmesh(nullptr),
  fe_coll(nullptr),
  scalar_fespace(nullptr),
  vector_fespace(nullptr),
  fespace(),
  u_0_gf()
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
                                    
                                    

        )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 2.0.0" << endl << LINE << endl;
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
    // Print NewtonSolverSettings
    if (user_input.UsingNewton())
        PrintNewtonSolverSettings();
    
    //----------------------------------------------------------------------
    // Print Output settings (if unsteady)
    if (user_input.UsingTimeIntegration())
        PrintOutput();

    //----------------------------------------------------------------------
    // Initialize material properties w/ solution field (only pertinent for non-uniform non-constant ones)
    if (rank == 0)
        cout << "\nInitializing material properties field...";

    UpdateMatProps(false);

    if (rank == 0)
        cout << " Done!" << endl;

    //---------------------------------------------------------------------
    // Create JOTSIterator object/s
    if (rank == 0)
    {
        cout << LINE << endl;
        cout << "Initializing iterator/s... ";
    }

    switch (Simulation_Type_Map.at(user_input.GetSimTypeLabel()))
    {
        case SIMULATION_TYPE::LINEARIZED_UNSTEADY:
            jots_iterator[PHYSICS_TYPE::THERMAL] = new LinearConductionOperator(user_input,
                                             boundary_conditions[PHYSICS_TYPE::THERMAL],
                                             all_bdr_attr_markers,
                                             *mat_props[MATERIAL_PROPERTY::DENSITY],
                                             *mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT],
                                             *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY],
                                             *fespace[PHYSICS_TYPE::THERMAL],
                                             time,
                                             dt);
            break;
        case SIMULATION_TYPE::NONLINEAR_UNSTEADY:
            jots_iterator[PHYSICS_TYPE::THERMAL] = new NonlinearConductionOperator(user_input,
                                             boundary_conditions[PHYSICS_TYPE::THERMAL],
                                             all_bdr_attr_markers,
                                             *mat_props[MATERIAL_PROPERTY::DENSITY],
                                             *mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT],
                                             *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY],
                                             *fespace[PHYSICS_TYPE::THERMAL],
                                             time,
                                             dt);
            break;
        case SIMULATION_TYPE::STEADY:
            jots_iterator[PHYSICS_TYPE::THERMAL] = new SteadyConductionOperator(user_input, 
                                                    boundary_conditions[PHYSICS_TYPE::THERMAL],
                                                    all_bdr_attr_markers,
                                                    *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY],
                                                    *fespace[PHYSICS_TYPE::THERMAL]);
            break;
    }
    //----------------------------------------------------------------------
    if (rank == 0)
        cout << "Done!" << endl;
}

void JOTSDriver::ProcessFiniteElementSetup()
{
    // Print appropriate solver type + get required material properties, add them to map
    if (rank == 0)
        cout << "Simulation Type: " << user_input.GetSimTypeLabel() << endl;

    // If not restart, refine mesh and initialize; else load VisItDataCollection
    if (!user_input.UsingRestart())
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
        // Define H1-conforming FE collection
        fe_coll = new H1_FECollection(user_input.GetFEOrder(), dim);
        
        // Initialize a scalar and vector ParFESpace using H1-conforming elements
        scalar_fespace = new ParFiniteElementSpace(pmesh, fe_coll, 1);
        vector_fespace = new ParFiniteElementSpace(pmesh, fe_coll, dim);

        // Initialize solution vector/s
        vector<string> soln_inits = user_input.GetInitializations();
        for (size_t i = 0; i < soln_inits.size(); i++)
        {   
            // Get associated physics type + solution name for each initialization
            PHYSICS_TYPE phys = Physics_Type_Map.at(soln_inits[i]);
            string solution_name = Solution_Names_Map.at(phys);

            // Initialize fespaces for this solution name / physics type
            switch (phys)
            {
                case PHYSICS_TYPE::THERMAL:
                    fespace[phys] = scalar_fespace;
                    break;
                case PHYSICS_TYPE::STRUCTURAL:
                    fespace[phys] = vector_fespace;
                    break;
                default:
                    // This would not be reached as code would crash at Physics_Type_Map.at() call
                    break;
            }
            
            u_0_gf[phys] = new ParGridFunction(fespace[phys]);
            ConstantCoefficient coeff_IC(user_input.GetInitialValue(soln_inits[i]));
            u_0_gf[phys]->ProjectCoefficient(coeff_IC);

            if (rank == 0)
                cout << "Initial " << solution_name << ": " << user_input.GetInitialValue(soln_inits[i]) << endl;
        }
        

    }
    //----------------------------------------------------------------------
    // Else using restart:
    else
    {   
        if (rank == 0)
        {
            cout << "Restart simulation..." << endl;
            cout << "Input Restart Root: " << user_input.GetRestartPrefix() << "_" << user_input.GetRestartCycle() << endl;
        }
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
        // By default, just check all solution names in DC
        for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
        {
            PHYSICS_TYPE phys = PHYSICS_TYPE(i);
            string soln_name = Solution_Names_Map.at(phys);
            u_0_gf[phys] = temp_visit_dc.GetParField(soln_name);
            
            // Skip if no PGF found
            if (!u_0_gf[phys])
                continue;

            if (!fe_coll)// Same FECollection used for thermal + structural -- H1-conforming
                fe_coll = u_0_gf[phys]->OwnFEC();

            // Save FESpace based on what type of physics it is modeling
            switch (phys)
            {
                case PHYSICS_TYPE::THERMAL:
                    scalar_fespace = u_0_gf[phys]->ParFESpace();
                    fespace[phys] = scalar_fespace;
                    break;
                case PHYSICS_TYPE::STRUCTURAL:
                    vector_fespace = u_0_gf[phys]->ParFESpace();
                    fespace[phys] = vector_fespace;
                    break;
                default:
                    break;
            }
            u_0_gf[phys]->MakeOwner(NULL);
            if (rank == 0)
                cout << soln_name << " Loaded!" << endl;
        }

        //  Ensure existence of both scalar and vector fespace
        if (!scalar_fespace)
            scalar_fespace = new ParFiniteElementSpace(pmesh, fe_coll, 1);
        if (!vector_fespace)
            vector_fespace = new ParFiniteElementSpace(pmesh, fe_coll, dim);

        // Now take ownership of mesh and fields
        // NOTE: If any further fields are added to restarts, they MUST be taken and deallocated here!!!!!
        // pmesh deallocated in destructor
        temp_visit_dc.SetOwnData(false);

        // Set the it_num and time
        it_num = temp_visit_dc.GetCycle();
        time = temp_visit_dc.GetTime();


        if (rank == 0)
        {
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
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
    {
        if (!fespace[i])
            continue;
        
        HYPRE_BigInt fe_size = fespace[i]->GlobalTrueVSize();
        if (rank == 0)
            cout << "Number of " << Solution_Names_Map.at(PHYSICS_TYPE(i)) << " DOFs: " << fe_size << endl;
    }
    //----------------------------------------------------------------------
    // Instantiate OutputManager
    output = new OutputManager(rank, user_input, *scalar_fespace);

    //----------------------------------------------------------------------
    // Set solution vectors from IC
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
    {
        if (!u_0_gf[i])
            continue;

        u[i] = new Vector();
        u_0_gf[i]->GetTrueDofs(*u[i]);
        // Add solution vector to OutputManager
        output->RegisterSolutionVector(Solution_Names_Map.at(PHYSICS_TYPE(i)), *u[i], *fespace[i]);
    }
    

}

void JOTSDriver::ProcessMaterialProperties()
{   
    if (rank == 0)
        cout << "\n";

    // Loop over INPUT material property info map from Config
    // Get keys
    vector<string> in_mat_prop_labels = user_input.GetMaterialPropertyLabels();

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

                mat_props[mp] = new PolynomialProperty(poly_coeffs, *scalar_fespace);
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
        cout << "Time Print Frequency: " << user_input.GetTimePrintFreq() << endl;
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

    precice_interface = new JOTSSolverInterface(user_input.GetPreciceParticipantName(), user_input.GetPreciceConfigFile(), rank, size, comm);

    if (precice_interface->getDimensions() != dim)
    {
        MFEM_ABORT("preCICE dimensions and mesh file dimensions are not the same!");
        return;
    }
}

void JOTSDriver::ProcessBoundaryConditions()
{
    // Verify input config BCs appropriately match input mesh BCs
    // Loop through user input BCs
    vector<string> in_bc_types = user_input.GetBCTypes();
    for (size_t i = 0; i < in_bc_types.size(); i++)
    {
        string type = in_bc_types[i];

        if (rank == 0)
            cout << "\n" << type << " Boundary Conditions:" << endl;

        vector<int> bc_keys = user_input.GetBCAttributes(type);

        // Confirm user input matches mesh bdr_attributes size...
        // Check count
        if (pmesh->bdr_attributes.Size() != bc_keys.size())
        {
            MFEM_ABORT("Input file BC count and mesh BC counts do not match.");
            return;
        }

        // Check one-to-oneness
        for (int j = 0; j < pmesh->bdr_attributes.Size(); j++)
        {
            int attr = pmesh->bdr_attributes[j];
            bool one_to_one = false;

            for (size_t k = 0; k < bc_keys.size(); k++)
            {
                if (attr == bc_keys[k])
                {
                    one_to_one = true;
                    k = bc_keys.size();
                }
            }
            if (!one_to_one)
            {   
                stringstream sstm;
                sstm << "No matching boundary attribute in config file for attribute " << attr;
                MFEM_ABORT(sstm.str());
                return;
            }
        }
        //----------------------------------------------------------------------
        // At this point, have verified exact number + one-to-one matching of BCs in config + mesh
        // (all bc_keys == pmesh->bdr_attributes essentially)
        // Instantiate boundary conditions + save any precice bc indices
        vector<int> precice_bc_indices;

        PHYSICS_TYPE phys = Physics_Type_Map.at(type);
        boundary_conditions[phys] = new BoundaryCondition*[pmesh->bdr_attributes.Size()];
        for (int j = 0; j < pmesh->bdr_attributes.Size(); j++)
        {   
            int attr = pmesh->bdr_attributes[j];
            vector<string> bc_info = user_input.GetBCInfo(type, attr);

            // If these are thermal boundary conditions
            if (phys == PHYSICS_TYPE::THERMAL)
            {
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
                        boundary_conditions[phys][j] =  new UniformConstantIsothermalBC(attr, value);
                        break;
                    case BOUNDARY_CONDITION::HEATFLUX:
                        value = stod(bc_info[1].c_str());
                        boundary_conditions[phys][j] =  new UniformConstantHeatFluxBC(attr, value);
                        break;
                    case BOUNDARY_CONDITION::PRECICE_ISOTHERMAL:
                        mesh_name = bc_info[1];
                        value = stod(bc_info[2].c_str());
                        boundary_conditions[phys][j] =  new PreciceIsothermalBC(attr, *fespace[PHYSICS_TYPE::THERMAL], mesh_name, value, *mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY]);
                        precice_bc_indices.push_back(j);
                        break;
                    case BOUNDARY_CONDITION::PRECICE_HEATFLUX:
                        mesh_name = bc_info[1];
                        value = stod(bc_info[2].c_str());
                        boundary_conditions[phys][j] =  new PreciceHeatFluxBC(attr, *fespace[PHYSICS_TYPE::THERMAL], mesh_name, value);
                        precice_bc_indices.push_back(j);
                        break;
                    case BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL:
                        amp = stod(bc_info[1].c_str());
                        afreq = stod(bc_info[2].c_str());
                        phase = stod(bc_info[3].c_str());
                        shift = stod(bc_info[4].c_str());
                        boundary_conditions[phys][j] = new UniformSinusoidalIsothermalBC(attr, amp, afreq, phase, shift);
                        break;
                    case BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX:
                        amp = stod(bc_info[1].c_str());
                        afreq = stod(bc_info[2].c_str());
                        phase = stod(bc_info[3].c_str());
                        shift = stod(bc_info[4].c_str());
                        boundary_conditions[phys][j] = new UniformSinusoidalHeatFluxBC(attr, amp, afreq, phase, shift);
                        break;
                    default:
                        MFEM_ABORT("Invalid/Unknown boundary condition specified: '" + bc_info[0] + "'");
                        return;
                }
            }
            else if (phys == PHYSICS_TYPE::STRUCTURAL)
            {
                MFEM_ABORT("Structural BCs not yet implemented");
                return;
            }
            else
            {
                MFEM_ABORT("Invalid/Unknown boundary condition type specified: '" + type + "'");
                return;
            }
        }
        //----------------------------------------------------------------------
        // Send precice bcs to precice_interface
        if (user_input.UsingPrecice())
            precice_interface->SetPreciceBCs(boundary_conditions[phys], precice_bc_indices);
        
        //----------------------------------------------------------------------
        // Print BCs
        for (int i = 0; i < pmesh->bdr_attributes.Size(); i++)
        {   
            BoundaryCondition* bc = boundary_conditions[phys][i];

            if (rank == 0)
            {
                cout << "Boundary Attribute " << bc->GetBdrAttr() << ": " << bc->GetInitString() << endl;
            }
        }
    }

    // Prepare BC attr arrays for applying coefficients
    all_bdr_attr_markers = new Array<int>[pmesh->bdr_attributes.Size()];
    for (int j = 0; j < pmesh->bdr_attributes.Size(); j++)
    {
        Array<int> bdr_attr(pmesh->bdr_attributes.Size());
        bdr_attr = 0;
        bdr_attr[j] = 1;
        all_bdr_attr_markers[j] = bdr_attr;
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
        cout << "Print Level: ";
        vector<string> pl = user_input.GetLinSolPrintLevel();
        for (size_t i = 0; i < pl.size(); i++)
        {
            cout << pl[i];
            if (i + 1 < pl.size())
                cout << ',';
            else
                cout << endl;
        }
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
        cout << "Print Level: ";
        vector<string> pl = user_input.GetNewtonPrintLevel();
        for (size_t i = 0; i < pl.size(); i++)
        {
            cout << pl[i];
            if (i + 1 < pl.size())
                cout << ',';
            else
                cout << endl;
        }

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
    if (user_input.UsingPrecice())
    {
        precice_dt = precice_interface->initialize();
        if (precice_interface->isActionRequired(JOTSSolverInterface::cowid))
        {
            // TODO: Only implemented for THERMAL problems presently
            u_0_gf[PHYSICS_TYPE::THERMAL]->SetFromTrueDofs(*u[PHYSICS_TYPE::THERMAL]);
            precice_interface->WriteData(*u_0_gf[PHYSICS_TYPE::THERMAL]);
            precice_interface->markActionFulfilled(JOTSSolverInterface::cowid);
        }
        precice_interface->initializeData();
    }
    
    while ((!user_input.UsingPrecice() && it_num < max_timesteps) 
    || (user_input.UsingPrecice() && precice_interface->isCouplingOngoing())) // use short-circuiting
    {
        // Handle preCICE calls
        if (user_input.UsingPrecice())
        {
            // Implicit coupling: save state
            if (precice_interface->isActionRequired(JOTSSolverInterface::cowic))
            {
                precice_saved_time = time;
                precice_saved_it_num = it_num;
                precice_interface->SaveOldState(*u[PHYSICS_TYPE::THERMAL]);
                precice_interface->markActionFulfilled(JOTSSolverInterface::cowic);
            }

            // Get read data
            if (precice_interface->isReadDataAvailable())
                precice_interface->GetReadData();

            // Update timestep if needed
            if (precice_dt < dt)
                dt = precice_dt;
        }

        // Update and Apply BCs
        UpdateAndApplyBCs();

        // Update MPs + apply changes through to iterator 
        // (note that some may have been changed by non-constant DBCs)
        UpdateMatProps(true);

        // Iterate
        Iteration();

        // Write any preCICE data, reload state if needed
        if (user_input.UsingPrecice())
        {
            if (precice_interface->isWriteDataRequired(dt))
            {
                u_0_gf[PHYSICS_TYPE::THERMAL]->SetFromTrueDofs(*u[PHYSICS_TYPE::THERMAL]);
                // Note: HF is calculated using k.GetLocalValue(), so no UpdateMatProps needed
                precice_interface->WriteData(*u_0_gf[PHYSICS_TYPE::THERMAL]);
            }
            // Advance preCICE
            precice_interface->advance(dt);


            // Implicit coupling
            if (precice_interface->isActionRequired(JOTSSolverInterface::coric))
            {
                time = precice_saved_time;
                it_num = precice_saved_it_num;
                precice_interface->ReloadOldState(*u[PHYSICS_TYPE::THERMAL]);
                precice_interface->markActionFulfilled(JOTSSolverInterface::coric);
                continue; // skip printing of timestep info AND outputting
            }
        }
   
        PostprocessIteration();
    }
    
    // Finalize preCICE if used
    if (user_input.UsingPrecice())
        precice_interface->finalize();


}

void JOTSDriver::UpdateMatProps(const bool apply_changes)
{
    for (int mp = 0; mp < MATERIAL_PROPERTY_SIZE; mp++)
    {
        if (mat_props[mp] != nullptr && !mat_props[mp]->IsConstant())
        {
            // Update coefficients using current solution field, ie: k=k(T)
            // TODO: just update using scalar thermal field for now
            mat_props[mp]->UpdateAllCoeffs(*u[PHYSICS_TYPE::THERMAL]);
            
            // Update any BLFs affected by changed coefficient (Apply)
            if (apply_changes)
            {
                for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
                {
                    if (!jots_iterator[i])
                        continue;
                    jots_iterator[i]->ProcessMatPropUpdate(MATERIAL_PROPERTY(mp));
                }
            }
        }
    }
}

void JOTSDriver::UpdateAndApplyBCs()
{
    // Update time-dependent/non-constant boundary condition coefficients
    // Get GF from current solution

    // Loop over every physics solver
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
    {
        if (!boundary_conditions[i])
            continue;
        
        bool n_changed = false;
        bool d_changed = false;

        u_0_gf[i]->SetFromTrueDofs(*u[i]);

        // Update BCs for this given physics solver
        for (int j = 0; j < pmesh->bdr_attributes.Size(); j++)
        {   
            // Update coefficients (could be preCICE calls, could be SetTime calls, etc.)
            if (!boundary_conditions[i][j]->IsConstant() || !initialized_bcs) // If not constant in time or not yet initialized for first iteration
            {
                boundary_conditions[i][j]->UpdateCoeff(time);

                if (boundary_conditions[i][j]->IsEssential())
                {
                    // Project correct values on boundary for essential BCs
                    d_changed = true;
                    u_0_gf[i]->ProjectBdrCoefficient(boundary_conditions[i][j]->GetCoeffRef(), all_bdr_attr_markers[j]);
                }
                else
                {
                    n_changed = true; // Flag that Neumann LF must be updated
                }
            }
        }

        if (d_changed)
        {
            // Apply changed Dirichlet BCs to u
            u_0_gf[i]->GetTrueDofs(*u[i]);
            // Here is where can set tmp_du_dt to be used if HO wanted
            // For higher-order, need T from n-1 timestep in addition to n timestep? is this worth doing?
            // Need to save previous timestep temperature in restarts
        }

        // If any non-constant Neumann terms, update Neumann linear form
        if (n_changed)
            jots_iterator[i]->UpdateNeumann();
    }
    if (!initialized_bcs)
        initialized_bcs = true;
}


void JOTSDriver::Iteration()
{
    // Step simulation
    // NOTE: Do NOT use ANY ODE-Solvers that update dt
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
        if (jots_iterator[i])
            jots_iterator[i]->Iterate(*u[i]);
    it_num++; // increment it_num - universal metric
}

void JOTSDriver::PostprocessIteration()
{   

    // Print timestep information if using TimeIntegration:
    bool time_out = user_input.UsingTimeIntegration() && it_num % user_input.GetTimePrintFreq() == 0;
    if (time_out)
    {
        if (rank == 0)
            printf("Step #%10i || Time: %1.6e out of %-1.6e || dt: %1.6e \n", it_num, time, dt*max_timesteps, dt);
    }
    
    // Check if blow up
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
    {
        if (!u[i])
            continue;
        if (u[i]->Max() > 1e10)
        {
            MFEM_ABORT("JOTS has blown up");
            return;
        }
    }

    // Write any output files if required
    // If unsteady simulation, then output at input frequencies
    // If steady, then output ONLY IF viz_freq and restart_freq != 0
    if (!user_input.UsingTimeIntegration())
    {
        // End of simulation reached. Update MPs prior to output w/ most recent solution field.
        UpdateMatProps(false);
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
            // Update MPs prior to output w/ most recent solution field.
            UpdateMatProps(false);
            if (rank == 0)
            {    
                cout << LINE << endl;
                cout << "Cycle: " << it_num << endl;
                cout << "Time: " << time << endl;
            }    
            if (viz_out)
            {
                if (rank == 0)
                    cout << "Saving Paraview Data..." <<endl;
                output->WriteVizOutput(it_num, time);
            }

            if (res_out)
            {
                if (rank == 0)
                    cout << "Saving Restart File..." << endl;
                output->WriteRestartOutput(it_num, time);
            }

            if (rank == 0)
                cout << LINE << endl;
        }
    }
}

JOTSDriver::~JOTSDriver()
{
    for (int i = 0; i < PHYSICS_TYPE_SIZE; i++)
    {   
        delete jots_iterator[i];
        delete u_0_gf[i];
        delete u[i];
        if (boundary_conditions[i])
        {
            for (int j = 0; j < pmesh->bdr_attributes.Size(); j++)
                delete boundary_conditions[i][j];
            delete[] boundary_conditions[i];
        }
        else
            delete boundary_conditions[i];
    }
    for (int i = 0; i < MATERIAL_PROPERTY_SIZE; i++)
        delete mat_props[i];
    delete precice_interface;
    delete[] all_bdr_attr_markers;
    delete pmesh;
    delete fe_coll;
    delete scalar_fespace;
    delete vector_fespace;
    delete output;

}