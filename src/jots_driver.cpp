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
                                
                                

    )" << endl << "MFEM-Based Thermal Solver w/ preCICE" << endl << "Version 1.0" << endl << line << endl;
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
    // Print the material properties and conductivity model
    if (myid == 0)
    {   
        cout << "\n\n";
        cout << "Density: " << user_input->GetDensity() << endl;
        cout << "Specific Heat Cp: " << user_input->GetCp() << endl;
        cout << "Thermal Diffusivity Model: ";
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
    double t_0;
    if (!user_input->UsesRestart()) // If not using restart
    {  
        //Create constant coefficient of initial temperature
        ConstantCoefficient T_0(user_input->GetInitialTemp());

        // Project that coefficient onto the GridFunction
        T_gf->ProjectCoefficient(T_0);

        t_0 = 0.0;
        if (myid == 0)
            cout << "Non-restart simulation --> Initial temperature field: " << user_input->GetInitialTemp() << endl;

    }
    else //else ARE using restart file
    {
        //TODO:
        if (myid == 0)
            cout << "Restarted simulation --> Restart file: " << user_input->GetRestartFile() << endl;

    }
    //----------------------------------------------------------------------
    // Print BCs - they will just be sent to ConductionOperator
    for (size_t i = 0; i < user_input->GetBCCount(); i++)
    {   
        BoundaryCondition* bc = user_input->GetBCs()[i];

        if (myid == 0)
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
    //----------------------------------------------------------------------
    // Create vector for holding true DOFs + instantiate ConductionOperator, sending all necessary parameters
    Vector u;
    T_gf->GetTrueDofs(u);
    oper = new ConductionOperator(user_input, *fespace, u, t_0);
}

void JOTSDriver::Run()
{
    
    //Here is where the run shit is done
    // preCICE may be implemented here, but boundary conditions are not fixed is an issue
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