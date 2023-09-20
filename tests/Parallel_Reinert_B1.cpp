#include <cmath>

#include "jots_driver.hpp"

#include "test_info.hpp"

using namespace std;
using namespace mfem;

double Reinert_B1_Analytical(const Vector& x, double time);

const int SIM_TIME = 2.0;
const double EPSILON = 1e-5;
const int N = 100; // Inclusive end of summation over terms of analytical solution infinite series to include

int main(int argc, char *argv[])
{  
    
    // Initialize MPI and HYPRE.
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // Get input file for Reinert_B1
    stringstream input_file;
    input_file << SOURCE_DIR << "/examples/Reinert_B1/Reinert_B1.ini";

    // Parse the config file
    Config input(input_file.str().c_str());

    // Update the mesh file location (may not be relative to that local directory)
    stringstream mesh_file;
    mesh_file << SOURCE_DIR << "/examples/meshes/Reinert_1D_block.mesh";
    input.SetMeshFile(mesh_file.str().c_str());

    // Suppress any output
    input.SetRestartFreq(INT_MAX);
    input.SetVisFreq(INT_MAX);

    // Only run for SIM_TIME
    input.SetFinalTime(SIM_TIME);

    // Reduce dt for test
    input.Setdt(1e-4);

    // Create new JOTSDriver
    JOTSDriver* driver = new JOTSDriver(input, myid, num_procs);

    // Run driver
    driver->Run();

    // Define exact solution coefficient
    FunctionCoefficient u_exact(Reinert_B1_Analytical);
    u_exact.SetTime(SIM_TIME);

    // Get the finite-element approximation solution
    double error = driver->GetOutputManager()->GetT_gf()->ComputeL2Error(u_exact);

    if (myid == 0)
        cout << "Error: " << error;

    if (error > EPSILON)
    {
        if (myid == 0)
            cout << " > " << EPSILON << endl << "Failed!" << endl;
        return 1;
    }
    else
        if (myid == 0)
            cout << " <= " << EPSILON << endl << "Success!" << endl;

    // Delete driver
    delete driver;
    

    return 0;
}

double Reinert_B1_Analytical(const Vector& x, double time)
{
    double theta = 1.0;
    double alpha = 2.5e-6;
    double L = 0.01;

    for (int n = 0; n < N+1; n++)
    {
        double A = pow(-1.0,n)/(2.0*n+1.0);
        double B = (-pow(2.0*n+1.0,2.0)*pow(M_PI,2)*alpha*time)/(4.0*pow(L,2));
        double C = ((2.0*n+1.0)*M_PI*x[0])/(2.0*L);

        theta = theta - (4.0/M_PI)*A*exp(B)*cos(C);
    }
    
    double T_0 = 300.0;
    double T_D = 500.0;

    return theta*(T_D - T_0) + T_0;

}