#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{
    const double EPSILON = 1e-9;

    // Initialize MPI and HYPRE.
    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // Get input file for Initial_Sinusoidal_Heated_Bar
    stringstream input_file_0;
    input_file_0 << SOURCE_DIR << "/examples/Restarted_Sinusoidal_Isothermal_Heated_Bar/Initial_Sinusoidal_Isothermal_Heated_Bar.ini";

    // Parse the coninput_file_0fig file
    Config input_0(input_file_0.str().c_str());

    // Update the mesh file location (may not be relative to that local directory)
    stringstream mesh_file;
    mesh_file << SOURCE_DIR << "/examples/meshes/bar.mesh";
    input_0.SetMeshFile(mesh_file.str().c_str());

    // Suppress any visualization output
    input_0.SetVisFreq(INT_MAX);

    // Only output at final iteration
    input_0.SetRestartFreq(input_0.GetFinalTime()/input_0.Getdt());

    // Create new JOTSDriver
    JOTSDriver* driver_0 = new JOTSDriver(input_0, myid, num_procs);

    // Run driver to end
    driver_0->Run();

    // Note that restart file should be outputted in whatever current WD is

    // Get input file for Restarted_Sinusoidal_Heated_Bar
    stringstream input_file_r;
    input_file_r << SOURCE_DIR << "/examples/Restarted_Sinusoidal_Isothermal_Heated_Bar/Restarted_Sinusoidal_Isothermal_Heated_Bar.ini";

    // Parse the config file
    Config input_r(input_file_r.str().c_str());

    // Create restarted JOTSDriver
    JOTSDriver* driver_r = new JOTSDriver(input_r, myid, num_procs);

    cout << "RESTART DRIVER CREATED" << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // Check difference between final temperature grid function of initial with restarted
    GridFunctionCoefficient u_0_final(driver_0->GetOutputManager()->GetT_gf());

    cout << "GFCOEFF CREATED" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // Get the finite-element approximation solution
    double error = driver_r->GetOutputManager()->GetT_gf()->ComputeL2Error(u_0_final);

    cout << "ERROR CALCULCATED" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0)
        cout << "Error: " << error;

    if (isnan(error))
    {
        cout << "Failed!" << endl;
        return 1;
    }

    if (error > EPSILON)
    {
        if (myid == 0)
            cout << " > " << EPSILON << endl << "Failed!" << endl;
        return 1;
    }
    else
        if (myid == 0)
            cout << " <= " << EPSILON << endl << "Success!" << endl;

    // Delete the drivers
    delete driver_0;
    delete driver_r;

    return 0;
}