#include "jots_driver.hpp"

using namespace std;

JOTSDriver::JOTSDriver(const char* input_file, int myid)
{   
    const std::string line = "-------------------------------------------------------------------------------------------";
        
    cout << line << endl;
    cout << R"(

         _  ____ _______ _____ 
        | |/ __ \__   __/ ____|
        | | |  | | | | | (___  
    _   | | |  | | | |  \___ \ 
   | |__| | |__| | | |  ____) |
    \____/ \____/  |_| |_____/ 
                                
                                

    )" << endl << line << endl;
    //----------------------------------------------------------------------
    // Parse config file
    user_input = new Config(input_file);
    if (myid == 0)
      cout << "Configuration file " << input_file << " parsed successfully!" << endl;
    //----------------------------------------------------------------------
    // Read the serial mesh
    const char* mesh_file = user_input->GetMeshFile().c_str();
    Mesh* mesh = new Mesh(mesh_file, 1);//pass one to generate edges
    dim = mesh->Dimension();
    if (myid == 0)
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
            if (myid == 0)
                cout << "Time Scheme: Euler Implicit" << endl;
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            ode_solver = new ForwardEulerSolver;
            if (myid == 0)
                cout << "Time Scheme: Euler Explicit" << endl;
            break;
        case TIME_SCHEME::RK4:
            ode_solver = new RK4Solver;
            if (myid == 0)
                cout << "Time Scheme: RK4" << endl;
            break;
        default:
            if (myid == 0)
                cout << "Unknown ODE solver type!" << endl;
            delete mesh;
            //TODO: throw error
    }
    //----------------------------------------------------------------------
    // Refine mesh in serial
    int ser_ref = user_input->GetSerialRefine();
    for (int lev = 0; lev < ser_ref; lev++)
    {
        mesh->UniformRefinement();
    }

    if (myid == 0)
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

    if (myid == 0)
        if (par_ref > 1)
            cout << "Parallel refinements of mesh completed: " << par_ref << endl;
        else
            cout << "No parallel refinements of mesh completed" << endl;

    //----------------------------------------------------------------------
    // Define parallel FE space on parallel mesh, and solution vector (grid function) T_gf
    if (user_input->GetFEOrder() > 1)
        fe_coll = new H1_FECollection(user_input->GetFEOrder(), dim);
    else
        fe_coll = pmesh->GetNodes()->OwnFEC();

    fespace = new ParFiniteElementSpace(pmesh, fe_coll);

    HYPRE_BigInt fe_size = fespace->GlobalTrueVSize();

    if (myid == 0)
    {
        cout << "Number of temperature unknowns: " << fe_size << endl;
    }

    T_gf = new ParGridFunction(fespace);
    //----------------------------------------------------------------------
    // Set the initial condition
    double t_0;
    if (!user_input->UsesRestart()) // If not using restart
    {  
        //Create constant coefficient of initial temperature
        ConstantCoefficient T_0(user_input->GetInitialTemp());

        // Project that coefficient onto the GridFunction
        T_gf->ProjectCoefficient(T_0);

        t_0 = 0.0;
        if (myid == 0)
            cout << "\n\nNon-restart simulation; initial temperature field: " << user_input->GetInitialTemp() << endl;

    }
    else //else ARE using restart file
    {
        //TODO:
        if (myid == 0)
            cout << "\n\nRestarted simulation. Not yet implemented! Closing..." << endl;
        //return 3;
    }
    //----------------------------------------------------------------------
    // Set up all BCs - note that these must be allowed to be dynamic
    //boundary_conditions = new BoundaryConditions(user_input->GetBCs());

    //----------------------------------------------------------------------
    // Instantiate ConductionOperator finally
    //oper = new ConductionOperator();
}

void JOTSDriver::Run()
{
    
    //Here is where the run shit is done
    // preCICE may be implemented here, but boundary conditions are not fixed is an issue
}