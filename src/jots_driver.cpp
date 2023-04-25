#include "jots_driver.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

JOTSDriver::JOTSDriver(const char* input_file, int myid, int num_procs)
{   
    rank = myid;
    size = num_procs;
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
      cout << "Configuration file: " << input_file << endl;
    //----------------------------------------------------------------------
    user_input->ReadFESetup();
    // Read the serial mesh
    const char* mesh_file = user_input->GetMeshFile().c_str();
    Mesh* mesh = new Mesh(mesh_file, 1);//pass one to generate edges
    dim = mesh->Dimension();
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

    if (rank == 0)
        cout << "Serial refinements of mesh completed: " << ser_ref << endl;


    //----------------------------------------------------------------------
    // Complete parallel decomposition of serial mesh, delete it, then refine parallel mesh
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    int par_ref = user_input->GetParallelRefine();

    for (int lev = 0; lev < par_ref; lev++)
    {
        pmesh->UniformRefinement();
    }

    if (rank == 0)
        cout << "Parallel refinements of mesh completed: " << par_ref << endl;


    if (rank == 0)
    {
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
    // Define parallel FE space on parallel mesh, and solution vector (grid function) T_gf
    fe_coll = new H1_FECollection(user_input->GetFEOrder(), dim);

    fespace = new ParFiniteElementSpace(pmesh, fe_coll);
    
    HYPRE_BigInt fe_size = fespace->GlobalTrueVSize();

    if (rank == 0)
    {
        cout << "Number of temperature nodes: " << fe_size << endl;
    }

    T_gf = new ParGridFunction(fespace);
    //----------------------------------------------------------------------
    user_input->ReadAndInitMatProps();
    // Print the material properties and conductivity model
    if (rank == 0)
    {   
        cout << "\n";
        cout << "Density: " << user_input->GetDensity() << endl;
        cout << "Specific Heat Cp: " << user_input->GetCp() << endl;
        cout << "Thermal Conductivity Model: " << user_input->GetConductivityModel()->GetInitString() << endl;
    }
    //----------------------------------------------------------------------
    user_input->ReadIC();
    if (rank == 0)
        cout << "\n";
    // Set the initial condition
    if (!user_input->UsesRestart()) // If not using restart
    {  
        //Create constant coefficient of initial temperature
        ConstantCoefficient T_0(user_input->GetInitialTemp());

        // Project that coefficient onto the GridFunction
        T_gf->ProjectCoefficient(T_0);

        if (rank == 0)
            cout << "Non-restart simulation --> Initial temperature field: " << user_input->GetInitialTemp() << endl;

    }
    else //else ARE using restart file
    {
        //TODO:
        if (rank == 0)
            cout << "Restarted simulation --> Restart file: " << user_input->GetRestartFile() << endl;

    }
    //----------------------------------------------------------------------
    // Read preCICE info
    user_input->ReadpreCICE();

    // Instantiate SolverInterface if using preCICE
    if (user_input->UsingpreCICE())
    {
        if (rank == 0)
        {
            cout << "\n";
            cout << "Using preCICE!" << endl;
            cout << "preCICE Participant Name: " << user_input->GetpreCICEParticipantName() << endl;
            cout << "preCICE Config File: " << user_input->GetpreCICEConfigFile() << endl;
        }
        interface = new SolverInterface(user_input->GetpreCICEParticipantName(), user_input->GetpreCICEConfigFile(), rank, size);

        if (interface->getDimensions() != dim)
            MFEM_ABORT("preCICE dimensions and mesh file dimensions are not the same!");
    
    }
    else
        interface = nullptr;
    //----------------------------------------------------------------------
    // Setup BCs
    if (user_input->UsingpreCICE())
        user_input->ReadAndInitBCs(T_gf, interface);
    else
        user_input->ReadAndInitBCs();

    if (rank == 0)
        cout << "\n";
    // Confirm user input matches mesh bdr_attributes...
    // Check count
    if (pmesh->bdr_attributes.Size() != user_input->GetBCCount())
    {
        if (rank == 0)
            cout << "Error: Input file BC count and mesh file BC counts do not match." << endl;
        return; // TODO: Error handling
    }
    // Check one-to-oneness
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {
        int attr = pmesh->bdr_attributes[i];
        bool one_to_one = true;

        for (size_t j = 0; j < user_input->GetBCCount(); j++)
        {
            if (attr == user_input->GetBCs()[j]->GetBdrAttr())
            {
                one_to_one = true;
                j = user_input->GetBCCount();
            }
        }
        if (!one_to_one)
        {   
            if (rank == 0)
                cout << "Error: No matching boundary attribute in mesh file for attribute " << attr << endl;
            return;// TODO: Error handling
        }
    }

    // Reorder BC array in Config to match bdr_attributes (ensures consistency when setting them)
    user_input->ReorderBCs(pmesh->bdr_attributes);
    // Print BCs - they will just be sent to ConductionOperator
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {   
        BoundaryCondition* bc = user_input->GetBCs()[i];

        if (rank == 0)
        {
            cout << "Boundary Attribute " << bc->GetBdrAttr() << ": " << bc->GetInitString() << endl;
        }
    }
    //----------------------------------------------------------------------
    user_input->ReadTimeInt();
    // Set ODE time integrator
    ode_solver = user_input->GetODESolver();
    if (rank == 0)
    {
        cout << "\n";
        cout << "Time Scheme: " << user_input->GetTimeSchemeString() << endl;
        cout << "Time Step: " << user_input->Getdt() << endl;
        cout << "Max Time: " << user_input->GetFinalTime() << endl;
    
    }
    //----------------------------------------------------------------------
    user_input->ReadLinSolSettings();
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
    user_input->ReadOutput();
    // Print output settings
    if (rank == 0)
    {
        cout << "\n\n";
        cout << "Restart Frequency: " << user_input->GetRestartFreq() << endl;
        cout << "Visualization Frequency: " << user_input->GetVisFreq() << endl;
    }
    //----------------------------------------------------------------------
    // Declare vector for holding true DOFs + instantiate ConductionOperator, sending all necessary parameters
    T_gf->GetTrueDofs(T);
    if (rank == 0)
    {
        cout << line << endl;
        cout << "Initializing operator... ";
    }
    oper = new ConductionOperator(user_input, *fespace, user_input->GetStartTime());// Needs BCs, FESpace, and initial time
    if (rank == 0)
        cout << "Done!" << endl;
}

void JOTSDriver::Run()
{   
    double time = user_input->GetStartTime();
    double dt = user_input->Getdt();
    double tf = user_input->GetFinalTime();
    int it_num = 0;


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

    double preCICE_dt = 0;

    // Initialize preCICE
    if (user_input->UsingpreCICE())
        preCICE_dt = interface->initialize();

    while (time < tf)//Main Solver Loop
    {
        // Apply the BCs + calculate thermal conductivities
        oper->PreprocessIteration(T, time);

        // Output IC:
        if (it_num == 0)
        {   
            T_gf->SetFromTrueDofs(T);
            paraview_dc.SetCycle(it_num);
            paraview_dc.SetTime(time);
            paraview_dc.RegisterField("Temperature",T_gf);
            paraview_dc.Save();
        }

        // Step in time - time automatically updated
        ode_solver->Step(T, time, dt);
        it_num += 1;

        // Print current timestep information:
        if (rank == 0)
            printf("Step #%10i || Time: %10.5g out of %-10.5g || dt: %10.5g \n", it_num, time, tf,  dt);
            //|| Rank 0 Max Temperature: %10.3g \n", it_num, time, tf,  dt, T.Max());
            //cout << "Step #" << it_num << " || t = " << time << "||" << "Rank 0 Max T: " << T.Max() << endl;
        
        if (T.Max() > 1e10)
        {
            MFEM_ABORT("JOTS has blown up");
            return;
        }

        if (it_num % user_input->GetVisFreq() == 0) // TODO: VisFreq must be nonzero
        {
            if (rank == 0)
                cout << line << endl << "Saving Paraview Data..." << endl << line << endl;
            // Save data in the ParaView format
            T_gf->SetFromTrueDofs(T);
            paraview_dc.SetCycle(it_num);
            paraview_dc.SetTime(time);
            paraview_dc.Save();
        }
    }

}

JOTSDriver::~JOTSDriver()
{
    delete interface;
    delete user_input;
    delete ode_solver;
    delete pmesh;
    delete fe_coll;
    delete T_gf;
    delete oper;

}