#include "jots_driver.hpp"

#include "../test_helper.hpp"

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
    string input_file_0 = "Restart_Files.ini";

    // Parse the config file
    Config input_0(input_file_0);

    // Create new JOTSDriver
    JOTSDriver* driver_0 = new JOTSDriver(input_0, myid, num_procs);

    // Run driver to end
    driver_0->Run();

    // Wait for all writings to end amongst all ranks
    MPI_Barrier(MPI_COMM_WORLD);

    // Note that restart file should be outputted in whatever current WD is

    // Get new Config and update to make it a restarted version of previous
    Config input_r(input_file_0);
    input_r.EnableRestart(true);
    input_r.SetRestartCycle(input_0.GetRestartFreq());
    
    // Create restarted JOTSDriver
    JOTSDriver* driver_r = new JOTSDriver(input_r, myid, num_procs);

    // Get the error per rank 
    // Not using ComputeL2Error because requires both GFs to be pointing to same Mesh in memory
    // ^(Not trivial)
    double error = driver_r->GetOutputManager()->GetVectorPGF("Temperature")->DistanceTo(*driver_0->GetOutputManager()->GetVectorPGF("Temperature"));

    // Print the error per rank
    for (int i = 0; i < num_procs; i++)
    {
        if (myid == i)
        {
            cout << "Rank " << i << " ";
            if (isnan(error))
            {   
                cout << "NaN detected -- Failed!" << endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            cout << "Error: " << error;

            if (error > EPSILON)
            {
                cout << " > " << EPSILON << endl << "Failed!" << endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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
    
    return 0;
}