#include "test_helper.hpp"

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Set rank and size
    int commRank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // Ensure that >1 proc being run for this script
    if (!(commSize > 1))
    {
        cout << "This test must be run with at least 2 MPI processes!" << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return 1;
    }

    // Split MPI Comms -- set Rank 0 to color 0, all others to color 1
    // Color 1 will run JOTS; Color 0 (Rank 0) will run dummy
    int color;
    if (commRank == 0)
    {
        color = 0;
    }
    else
    {
        color = 1;
    }
    MPI_Comm subComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, commRank, &subComm);

    // Get rank + size in new communicator
    int subRank = -1;
    MPI_Comm_rank(subComm, &subRank);

    int subSize = -1;
    MPI_Comm_size(subComm, &subSize);

    int exit;
    
    if (color == 0)
    {
        // Create Dummy Mesh as single point @ x = 0.0
        vector<double> vertices(3);
        vertices.at(0) = 0;
        vertices.at(1) = 0;
        vertices.at(2) = 0;

        // Create lambda function for constant send of 7.5e5 HF
        function<double(const Vector&, double)> Constant_HF = [=](const Vector&x, double t) -> double { return -7.5e5; };
        
        exit = preCICE_Dummy_Analytical_Test("Reinert_B3_preCICE_config.xml", "Temperature", "Heat-Flux", vertices, Constant_HF, Reinert_B3_Analytical, 100, subRank, subSize, subComm);
        if (exit == 1)
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else
    {
        // Run JOTS
        exit = RunTestConfig("Reinert_B3_preCICE.ini", subRank, subSize, subComm);
        if (exit == 1)
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Finalize MPI
    MPI_Finalize();

    // Return exit
    return exit;
}