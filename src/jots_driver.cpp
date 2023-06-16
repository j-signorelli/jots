#include "jots_driver.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

const string JOTSDriver::LINE = "-------------------------------------------------------------------";
const double JOTSDriver::TIME_TOLERANCE = 1e-12;

JOTSDriver::JOTSDriver(const char* input_file, const int myid, const int num_procs)
: rank(myid),
  size(num_procs),
  adapter(nullptr),
  user_input(nullptr),
  boundary_conditions(nullptr),
  cond_model(nullptr),
  ode_solver(nullptr),
  pmesh(nullptr),
  fe_coll(nullptr),
  fespace(nullptr),
  oper(nullptr)
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
                                    
                                    

        )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.0" << endl << LINE << endl;
    }
    //----------------------------------------------------------------------
    // Parse config file
    user_input = new Config(input_file);
    if (rank == 0)
      cout << "Configuration file " << input_file << " parsed successfully!" << endl;
    //----------------------------------------------------------------------
    // Create solution GF
    ParGridFunction* temp_T_gf;

    //----------------------------------------------------------------------
    // If not restart, refine mesh and initialize; else load VisItDataCollection
    if (!user_input->UsesRestart())
    {
        if (rank == 0)
            cout << "Non-restart simulation..." << endl;
        //----------------------------------------------------------------------
        // Create serial mesh
        const char* mesh_file = user_input->GetMeshFile().c_str();
        Mesh* mesh = new Mesh(mesh_file, 1);//pass one to generate edges
        dim = mesh->Dimension();
        //----------------------------------------------------------------------
        // Print mesh info
        if (rank == 0)
        {
            cout << "\n";
            cout << "Mesh File: " << mesh_file << endl;
        }
        //----------------------------------------------------------------------
        // Refine mesh in serial
        int ser_ref = user_input->GetSerialRefine();
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
        pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
        delete mesh;
        int par_ref = user_input->GetParallelRefine();

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
        fe_coll = new H1_FECollection(user_input->GetFEOrder(), dim);

        fespace = new ParFiniteElementSpace(pmesh, fe_coll);
        
        //----------------------------------------------------------------------
        // Create solution GF
        temp_T_gf = new ParGridFunction(fespace);
        
        //----------------------------------------------------------------------
        // Print initial condition info + set T_gf
        ConstantCoefficient coeff_IC(user_input->GetInitialTemp());
        temp_T_gf->ProjectCoefficient(coeff_IC);
        
        it_num = 0;
        time = 0.0;

        if (rank == 0)
            cout << "Initial temperature field: " << user_input->GetInitialTemp() << endl;

    }
    //----------------------------------------------------------------------
    // Else using restart:
    else
    {   
        if (rank == 0)
            cout << "Restart simulation..." << endl;
        //----------------------------------------------------------------------
        // Read in VisItDataCollection
        VisItDataCollection temp_visit_dc(MPI_COMM_WORLD, user_input->GetRestartPrefix());
        temp_visit_dc.Load(user_input->GetRestartCycle());

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
            cout << "Input Restart Root: " << user_input->GetRestartPrefix() << "_" << user_input->GetRestartCycle() << endl;
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
    // Print material properties
    if (rank == 0)
    {   
        cout << "\n";
        cout << "Density: " << user_input->GetDensity() << endl;
        cout << "Specific Heat Cp: " << user_input->GetCp() << endl;
    }
    //----------------------------------------------------------------------
    // Create ConductivityModel object
    vector<string> cond_info = user_input->GetCondInfo();
    vector<double> poly_coeffs;
    switch (Conductivity_Model_Map.at(cond_info[0]))
    {
        case CONDUCTIVITY_MODEL::UNIFORM: // Uniform conductivity
            cond_model = new UniformCond(stod(cond_info[1].c_str()));
            break;
        case CONDUCTIVITY_MODEL::POLYNOMIAL: // Polynomial model for conductivity

            for (int i = 1; i < cond_info.size(); i++)
                poly_coeffs.push_back(stod(cond_info[i].c_str()));

            cond_model = new PolynomialCond(poly_coeffs, *fespace);
            break;
        default:
            MFEM_ABORT("Unknown/Invalid thermal conductivity model specified");
            return;
    }
    //----------------------------------------------------------------------
    // Print conductivity model info
    if (rank == 0)
        cout << "Thermal Conductivity Model: " << cond_model->GetInitString() << endl;
    //----------------------------------------------------------------------
    // Set ODE time integrator
    ode_solver = user_input->GetODESolver();
    //----------------------------------------------------------------------
    // Print time integration information
    if (rank == 0)
    {
        cout << "\n";
        cout << "Time Scheme: " << user_input->GetTimeSchemeString() << endl;
        cout << "Time Step: " << user_input->Getdt() << endl;
        cout << "Max Time: " << user_input->GetFinalTime() << endl;
    
    }
    //----------------------------------------------------------------------
    // Set time stuff
    dt = user_input->Getdt();
    tf = user_input->GetFinalTime();
    //----------------------------------------------------------------------
    // Print any precice info + instantiate adapter object if needed
    if (user_input->UsingPrecice())
    {
        if (rank == 0)
        {
            cout << "\n";
            cout << "Using preCICE!" << endl;
            cout << "preCICE Participant Name: " << user_input->GetPreciceParticipantName() << endl;
            cout << "preCICE Config File: " << user_input->GetPreciceConfigFile() << endl;
        }

        adapter = new PreciceAdapter(user_input->GetPreciceParticipantName(), user_input->GetPreciceConfigFile(), rank, size);

        if (adapter->GetDimension() != dim)
        {
            MFEM_ABORT("preCICE dimensions and mesh file dimensions are not the same!");
            return;
        }
    }
    //----------------------------------------------------------------------
    // Verify input config BCs appropriately match input mesh BCs
    if (rank == 0)
        cout << "\n";

    // Confirm user input matches mesh bdr_attributes...
    // Check count
    if (pmesh->bdr_attributes.Size() != user_input->GetBCCount())
    {
        if (rank == 0)
            MFEM_ABORT("Input file BC count and mesh BC counts do not match.");
        return;
    }
    
    // Check one-to-oneness
    // Also create an array of indices that match mesh bdr attributes to input bdr attributes
    // ^ in event that user inputs bdr attributes in different order than mesh
    vector<int> bdr_index;
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {
        int attr = pmesh->bdr_attributes[i];
        bool one_to_one = false;

        for (size_t j = 0; j < user_input->GetBCCount(); j++)
        {   
            pair<int, vector<string>> bc = user_input->GetBCInfo(j);
            if (attr == bc.first)
            {
                one_to_one = true;
                bdr_index.push_back(j);
                j = user_input->GetBCCount();
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
    boundary_conditions = new BoundaryCondition*[user_input->GetBCCount()];
    for (int i = 0; i < user_input->GetBCCount(); i++)
    {
        pair<int, vector<string>> bc = user_input->GetBCInfo(bdr_index[i]);
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
    if (user_input->UsingPrecice())
        adapter->AddPreciceBCs(boundary_conditions, precice_bc_indices);

    //----------------------------------------------------------------------
    // Print BCs
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {   
        BoundaryCondition* bc = boundary_conditions[i];

        if (rank == 0)
        {
            cout << "Boundary Attribute " << bc->GetBdrAttr() << ": " << bc->GetInitString() << endl;
        }
    }
    //----------------------------------------------------------------------
    // Print linear solver settings
    if (rank == 0)
    {
        cout << "\n";
        cout << "Linear Solver: " << user_input->GetSolverString() << endl;
        cout << "Preconditioner: " << user_input->GetPrecString() << endl;
        cout << "Max Iterations: " << user_input->GetMaxIter() << endl;
        cout << "Absolute Tolerance: " << user_input->GetAbsTol() << endl;
        cout << "Relative Tolerance: " << user_input->GetRelTol() << endl;
    }
    //----------------------------------------------------------------------
    // Print output settings
    if (rank == 0)
    {
        cout << "\n\n";
        cout << "Restart Frequency: " << user_input->GetRestartFreq() << endl;
        cout << "Visualization Frequency: " << user_input->GetVisFreq() << endl;
    }
    //---------------------------------------------------------------------
    // Create main solution vector from IC
    temp_T_gf->GetTrueDofs(T);
    delete temp_T_gf;
    //----------------------------------------------------------------------
    // Instantiate ConductionOperator, sending all necessary parameters
    if (rank == 0)
    {
        cout << LINE << endl;
        cout << "Initializing operator... ";
    }
    oper = new ConductionOperator(user_input, boundary_conditions, cond_model, *fespace, time);
    if (rank == 0)
        cout << "Done!" << endl;
    //----------------------------------------------------------------------
    // Instantiate OutputManager
    output = new OutputManager(rank, fespace, user_input, T, cond_model);
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
    {
        output->WriteVizOutput(it_num, time);
    }


    // If using precice: initialize + get/send initial data if needed
    // Make actual precice interface calls directly here for readability/clarity
    if (user_input->UsingPrecice())
    {
        precice_dt = adapter->Interface()->initialize();
        if (adapter->Interface()->isActionRequired(PreciceAdapter::cowid))
        {
            adapter->WriteData(T, cond_model);
            adapter->Interface()->markActionFulfilled(PreciceAdapter::cowid);
        }
        adapter->Interface()->initializeData();
    }
    
    bool viz_out = false;
    bool res_out = false;

    while ( (!user_input->UsingPrecice() && time < tf) 
        || (user_input->UsingPrecice() && adapter->Interface()->isCouplingOngoing()))//Main Solver Loop - use short-circuiting
    {

        if (user_input->UsingPrecice())
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


        // Apply the BCs to state + any updated thermal conductivities
        oper->PreprocessIteration(T);


        // Update timestep if needed
        if (user_input->UsingPrecice())
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

        if (user_input->UsingPrecice())
        {
            if (adapter->Interface()->isWriteDataRequired(dt))
                adapter->WriteData(T, cond_model);

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
        viz_out = user_input->GetVisFreq() != 0 && it_num % user_input->GetVisFreq() == 0;
        res_out = user_input->GetRestartFreq() != 0 && it_num % user_input->GetRestartFreq() == 0;
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
    if (user_input->UsingPrecice())
        adapter->Interface()->finalize();


}

void JOTSDriver::UpdateMatProps()
{
    if (!cond_model->IsConstant())
        cond_model->UpdateCoeff(T); // Update k coefficient
}

void JOTSDriver::PreprocessIteration()
{
    // Update BCs

    // If conductivity not constant, it may have changed -> Update stiffness
    if (!cond_model->IsConstant())
        oper->UpdateStiffness();
}

JOTSDriver::~JOTSDriver()
{
    delete adapter;
    delete cond_model;
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
        delete boundary_conditions[i];
    delete[] boundary_conditions;
    delete user_input;
    delete ode_solver;
    delete pmesh;
    delete fe_coll;
    delete oper;
    delete output;

}