#include "jots_driver.hpp"

#include "test_info.hpp"

using namespace std;
using namespace mfem;

double Reinert_B1_Analytical(const Vector& x, double time);

const int SIM_TIME = 1.0;
const double EPSILON = 1e-9

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
    mesh_file << SOURCE_DIR << "/examples/meshes/block.mesh";
    input.SetMeshFile(mesh_file.str().c_str());

    // Suppress any output
    input.SetRestartFreq(INT_MAX);
    input.SetVisFreq(INT_MAX);

    // Only run for 1 physical second
    input.SetFinalTime(1.0);

    // Create new JOTSDriver
    JOTSDriver* driver = new JOTSDriver(input, myid, num_procs);

    // Run driver
    driver->Run();

    // Define exact solution coefficient
    FunctionCoefficient u_exact(Reinert_B1_Analytical);

    // Get the finite-element approximation solution
    double error = driver->GetOutputManager()->ComputeL2Error(u_exact);

    if (rank == 0)
        cout << "Error: " << error << endl;

    if (error > epsilon)
    {
        cout << "Failed!" << endl;
        return 1;
    }
    else
        cout << "Success!" << endl;

    // Delete driver
    delete driver;


    return 0;
}

double Reinert_B1_Analytical(const Vector& x, double time)
{

}