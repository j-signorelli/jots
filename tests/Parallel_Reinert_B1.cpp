#include "jots_driver.hpp"


using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{  
    // Initialize MPI and HYPRE.
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // Get input file for Reinert_B1
    const char *input_file = "";

    // Parse the config file
    Config input = new Config(input_file);


    // Create new JOTSDriver
    JOTSDriver* driver = new JOTSDriver(input, myid, num_procs);

    // Run driver
    driver->Run();

    // Delete driver
    delete driver;


    return 0;
}

