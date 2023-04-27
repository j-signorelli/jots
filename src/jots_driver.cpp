#include "jots_driver.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

JOTSDriver::JOTSDriver(const char* input_file, const int myid, const int num_procs)
: rank(myid),
  size(num_procs),
  adapter(nullptr),
  user_input(nullptr),
  boundary_conditions(nullptr),
  cond_model(nullptr),
  state(nullptr),
  ode_solver(nullptr),
  pmesh(nullptr),
  fe_coll(nullptr),
  fespace(nullptr),
  oper(nullptr)
{   
    if (rank == 0)
    {
        cout << line << endl;
        cout << R"(

              _  ____ _______ _____ 
             | |/ __ \__   __/ ____|
             | | |  | | | | | (___  
         _   | | |  | | | |  \___ \ 
        | |__| | |__| | | |  ____) |
         \____/ \____/  |_| |_____/ 
                                    
                                    

        )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.0" << endl << line << endl;
    }
    //----------------------------------------------------------------------
    // Parse config file
    user_input = new Config(input_file);
    if (rank == 0)
      cout << "Configuration file " << input_file << " parsed successfully!" << endl;
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
        cout << "Problem Dimension: " << dim << endl;
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
    // Print parallel refinements and mesh attributes from mesh file
    if (rank == 0)
    {
        cout << "Parallel refinements of mesh completed: " << par_ref << endl;
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
    // Define parallel FE space on parallel mesh
    fe_coll = new H1_FECollection(user_input->GetFEOrder(), dim);

    fespace = new ParFiniteElementSpace(pmesh, fe_coll);
    
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
    switch (Conductivity_Model_Map.at(bc_info[0]))
    {
        /* TODO: Linear
        case CONDUCTIVITYMODEL::L
            break;
        */
        case CONDUCTIVITY_MODEL::UNIFORM: // Uniform conductivity
            cond_model = new UniformCond(stod(bc_info[1].c_str()));
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
    // Print initial condition info + prepare initial condition gridfunction
    ParGridFunction temp_IC(fespace);
    double start_time = 0;
    int it_num = 0;
    if (rank == 0)
        cout << "\n";
    if (!user_input->UsesRestart()) // If not using restart
    {
        temp_IC.ProjectCoefficient(user_input->GetInitialTemp());
        
        if (rank == 0)
            cout << "Non-restart simulation --> Initial temperature field: " << user_input->GetInitialTemp() << endl;

    }
    else //else ARE using restart file
    {
        if (rank == 0)
            cout << "Restarted simulation --> Restart file: " << user_input->GetRestartFile() << endl;

        // TODO: read in restart file and set initial condition GF
    
    }
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
    // Create SolverState object for holding simulation info
    state = new SolverState(temp_IC, it_num, start_time, user_input->Getdt(), user_input->GetFinalTime());
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
            MFEM_ABORT("Input file BC count and mesh file BC counts do not match.");
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
                sstm << "No matching boundary attribute in mesh file for attribute " << attr;
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

        switch (Boundary_Condition_Map.at(bc.second[0]))
        {
            case BOUNDARY_CONDITION::ISOTHERMAL:
                double uniform_value = stod(bc.second[1].c_str())
                boundary_conditions[i] =  new UniformIsothermalBC(bc.first, uniform_value);
                break;
            case BOUNDARY_CONDITION::HEATFLUX:
                double uniform_value = stod(bc.second[1].c_str())
                boundary_conditions[i] =  new UniformHeatFluxBC(bc.first, uniform_value);
                break;
            case BOUNDARY_CONDITION::PRECICE_ISOTHERMAL:
                string mesh_name = bc.second[1];
                double initial_value = stod(bc.second[2].c_str())
                boundary_conditions[i] =  new PreciceIsothermalBC(bc.first, fespace, mesh_name, user_input->UsesRestart(), initial_value);
                precice_bc_indices.push_back(i);
                break;
            case BOUNDARY_CONDITION::PRECICE_HEATFLUX:
                string mesh_name = bc.second[1];
                double initial_value = stod(bc.second[2].c_str())
                boundary_conditions[i] =  new PreciceHeatFluxBC(bc.first, fespace, mesh_name, user_input->UsesRestart(), initial_value);
                precice_bc_indices.push_back(i);
                break;
            default:
                MFEM_ABORT("Invalid/Unknown boundary condition specified");
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
    //----------------------------------------------------------------------
    // Instantiate ConductionOperator, sending all necessary parameters
    if (rank == 0)
    {
        cout << line << endl;
        cout << "Initializing operator... ";
    }
    oper = new ConductionOperator(user_input, boundary_conditions, cond_model, *fespace, state->GetTime());
    if (rank == 0)
        cout << "Done!" << endl;
}

void JOTSDriver::Run()
{   

    // Initialize the ODE Solver
    if (rank == 0)
        cout << "Initializing solver...";

    ode_solver->Init(*oper);
    
    if (rank == 0)
        cout << " Done!" << endl;

    // Set up Paraview
    ParaViewDataCollection paraview_dc("Output", pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(user_input->GetFEOrder());
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);

    double precice_dt = 0;

    // If using precice: initialize + get/send initial data
    if (user_input->UsingPrecice())
    {
        precice_dt = adapter->Initialize();
        adapter->SendInitialData();
        adapter->InitializeData();
    }

    while ( (!user_input->UsingpreCICE() && state->GetTime() < state->GetFinalTime()) 
        || (user_input->UsingpreCICE() && interface->isCouplingOngoing()))//Main Solver Loop - use short-circuiting
    {
        // Apply the BCs to state + calculate thermal conductivities
        oper->PreprocessIteration(state->GetTRef());

        // Output IC:
        if (state->GetItNum() == 0)
        {   // TODO: Output class probably
            paraview_dc.SetCycle(state->GetItNum());
            paraview_dc.SetTime(state->GetTime());
            paraview_dc.RegisterField("Temperature",state->GetGF());
            paraview_dc.Save();
        }

        // Advance preCICE if using it
        if (user_input->UsingpreCICE())
        {
            interface->advance(preCICE_dt);
            if (preCICE_dt < state->Getdt())
                state->Setdt(preCICE_dt);
        }

        // Step in time - time automatically updated
        // NOTE: Do NOT use ANY ODE-Solvers that update dt
        ode_solver->Step(state->GetTRef(), state->GetTimeRef(), state->GetdtRef());
        state->SetItNum(state->GetItNum() + 1);

        // Print current timestep information:
        if (rank == 0)
            printf("Step #%10i || Time: %10.5g out of %-10.5g || dt: %10.5g \n", state->GetItNum(), state->GetTime(), state->GetFinalTime(),  state->Getdt());
            //|| Rank 0 Max Temperature: %10.3g \n", it_num, time, tf,  dt, T.Max());
            //cout << "Step #" << it_num << " || t = " << time << "||" << "Rank 0 Max T: " << T.Max() << endl;
        
        if (state->GetTRef().Max() > 1e10)
        {
            MFEM_ABORT("JOTS has blown up");
            return;
        }

        if (state->GetItNum() % user_input->GetVisFreq() == 0) // TODO: VisFreq must be nonzero
        {
            if (rank == 0)
                cout << line << endl << "Saving Paraview Data..." << endl << line << endl;
            // Save data in the ParaView format
            state->UpdateGF();
            paraview_dc.SetCycle(state->GetItNum());
            paraview_dc.SetTime(state->GetTime());
            paraview_dc.Save();
        }
    }
    if (user_input->UsingpreCICE())
        interface->finalize();


}

JOTSDriver::~JOTSDriver()
{
    delete adapter;
    delete cond_model;
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
        delete boundary_conditions[i];
    delete[] boundary_conditions;
    delete user_input;
    delete state;
    delete ode_solver;
    delete pmesh;
    delete fe_coll;
    delete oper;

}