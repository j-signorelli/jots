#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{
    const double EPSILON = 1e-12;

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

    // Wait for all writings to end amongst all ranks
    MPI_Barrier(MPI_COMM_WORLD);

    // Note that restart file should be outputted in whatever current WD is

    // Get input file for Restarted_Sinusoidal_Heated_Bar
    stringstream input_file_r;
    input_file_r << SOURCE_DIR << "/examples/Restarted_Sinusoidal_Isothermal_Heated_Bar/Restarted_Sinusoidal_Isothermal_Heated_Bar.ini";

    // Parse the config file
    Config input_r(input_file_r.str().c_str());

    // Create restarted JOTSDriver
    JOTSDriver* driver_r = new JOTSDriver(input_r, myid, num_procs);

    // Get the error per rank 
    // Not using ComputeL2Error because requires both GFs to be pointing to same Mesh in memory
    // ^(Not trivial)
    double error = driver_r->GetOutputManager()->GetT_gf()->DistanceTo(*driver_0->GetOutputManager()->GetT_gf());

    bool success = true;

    // Print the error per rank
    for (int i = 0; i < num_procs; i++)
    {
        if (myid == i)
        {
            cout << "Rank " << i << " ";
            if (isnan(error))
            {   
                cout << "NaN detected -- Failed!" << endl;
                success = false;
            }
            cout << "Error: " << error;

            if (error > EPSILON)
            {
                cout << " > " << EPSILON << endl << "Failed!" << endl;
                success = false;
            }
            else
                cout << " <= " << EPSILON << endl << "Success!" << endl;
            fflush(stdout);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Delete the drivers
    delete driver_0;
    delete driver_r;
    
    // Return appropriately if failed or not
    if (!success)
        return 1;

    return 0;
}