#include "jots_driver.hpp"

using namespace std;

JOTSDriver::JOTSDriver(const char* input_file, int myid)
{   
    rank = myid;
    
    cout << line << endl;
    cout << R"(

         _  ____ _______ _____ 
        | |/ __ \__   __/ ____|
        | | |  | | | | | (___  
    _   | | |  | | | |  \___ \ 
   | |__| | |__| | | |  ____) |
    \____/ \____/  |_| |_____/ 
                                
                                

    )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.0" << endl << line << endl;
    //----------------------------------------------------------------------
    // Parse config file
    user_input = new Config(input_file);
    if (rank == 0)
      cout << "Configuration file " << input_file << " parsed successfully!" << endl;
    //----------------------------------------------------------------------
    // Read the serial mesh
    const char* mesh_file = user_input->GetMeshFile().c_str();
    Mesh* mesh = new Mesh(mesh_file, 1);//pass one to generate edges
    dim = mesh->Dimension();
    if (rank == 0)
    {
        cout << "Mesh File: " << mesh_file << endl;
        cout << "Problem Dimension: " << dim << endl;
    }
    //----------------------------------------------------------------------
    // Set ODE time integrator
    switch (user_input->GetTimeScheme())
    {
        case TIME_SCHEME::EULER_IMPLICIT:
            ode_solver = new BackwardEulerSolver;
            if (rank == 0)
                cout << "Time Scheme: Euler Implicit" << endl;
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            ode_solver = new ForwardEulerSolver;
            if (rank == 0)
                cout << "Time Scheme: Euler Explicit" << endl;
            break;
        case TIME_SCHEME::RK4:
            ode_solver = new RK4Solver;
            if (rank == 0)
                cout << "Time Scheme: RK4" << endl;
            break;
    }
    //----------------------------------------------------------------------
    // Refine mesh in serial
    int ser_ref = user_input->GetSerialRefine();
    for (int lev = 0; lev < ser_ref; lev++)
    {
        mesh->UniformRefinement();
    }

    if (rank == 0)
        if (ser_ref > 1)
            cout << "Serial refinements of mesh completed: " << ser_ref << endl;
        else
            cout << "No serial refinements of mesh completed" << endl;

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
        if (par_ref > 1)
            cout << "Parallel refinements of mesh completed: " << par_ref << endl;
        else
            cout << "No parallel refinements of mesh completed" << endl;

    //----------------------------------------------------------------------
    // Define parallel FE space on parallel mesh, and solution vector (grid function) T_gf
    fe_coll = new H1_FECollection(user_input->GetFEOrder(), dim);

    fespace = new ParFiniteElementSpace(pmesh, fe_coll);

    HYPRE_BigInt fe_size = fespace->GlobalTrueVSize();

    if (rank == 0)
    {
        cout << "Number of temperature unknowns: " << fe_size << endl;
    }

    T_gf = new ParGridFunction(fespace);
    //----------------------------------------------------------------------
    // Print the material properties and conductivity model
    if (rank == 0)
    {   
        cout << "\n\n";
        cout << "Density: " << user_input->GetDensity() << endl;
        cout << "Specific Heat Cp: " << user_input->GetCp() << endl;
        cout << "Thermal Conductivity Model: ";
        switch (user_input->GetConductivityModel()->GetModel())
        {
            case CONDUCTIVITY_MODEL::CONSTANT:
                cout << "Constant -- k: " << ((ConstantCond*)user_input->GetConductivityModel())->Getk() << endl;
                break;
            case CONDUCTIVITY_MODEL::LINEARIZED:
                cout << "Linearized -- k: " << ((LinearizedCond*)user_input->GetConductivityModel())->Getk() << "; alpha: " << ((LinearizedCond*)user_input->GetConductivityModel())->Getalpha() << endl;
        }
    }
    //----------------------------------------------------------------------
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
    // Print BCs - they will just be sent to ConductionOperator
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {   
        BoundaryCondition* bc = user_input->GetBCs()[i];

        if (rank == 0)
        {
            cout << "Boundary Attribute " << i+1 << ": ";
            switch (bc->GetType())
            {
                case BOUNDARY_CONDITION::HEATFLUX:
                    cout << "Heat Flux";
                    break;
                case BOUNDARY_CONDITION::ISOTHERMAL:
                    cout << "Isothermal";
                    break;
            };

            cout << " --- Value: " << bc->GetValue() << endl;
        }
    }
    // Verify that listed BCs match the boundary attributes in the mesh
    
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
    double t0 = user_input->GetStartTime();
    double dt = user_input->Getdt();
    double tf = user_input->GetFinalTime();

    // Initialize the ODE Solver
    if (rank == 0)
        cout << "Initializing solver...";
    ode_solver->Init(*oper);
    
    if (rank == 0)
        cout << " Done!" << endl;

    double time = t0;
    int it_num = 0;

    // Set up Paraview
    ParaViewDataCollection paraview_dc("Output", pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(user_input->GetFEOrder());
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(t0);
    paraview_dc.RegisterField("Temperature",T_gf);
    paraview_dc.Save();

    while (time < tf)//Main Solver Loop- TODO: Fix this
    {

        // Apply the BCs
        oper->ApplyBCs(T);

        // Calculate thermal conductivities
        oper->SetThermalConductivities(T);

        // Step in time - time automatically updated
        ode_solver->Step(T, time, dt);
        it_num += 1;

        // Print current timestep information:
        if (rank == 0)
            printf("Step #%10i || Time: %10.5g out of %-10.5g || dt: %10.5g || Rank 0 Max Temperature: %10.3g \n", it_num, time, tf,  dt, T.Max());
            //cout << "Step #" << it_num << " || t = " << time << "||" << "Rank 0 Max T: " << T.Max() << endl;
        
	if (T.Max() > 1e10)
	{
        if (rank == 0)
	        cout << "Error: JOTS has diverged" << endl;
	    return; // TODO: Have main return nonzero value
	}

        if (it_num % user_input->GetVisFreq() == 0)
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
    delete user_input;
    delete ode_solver;
    delete pmesh;
    delete fe_coll;
    delete T_gf;
    delete oper;

}