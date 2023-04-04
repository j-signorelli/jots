//#include <fstream>
//#include <iostream>

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

   // Parse input file
   const char *input_file = "../config_template.ini";
   //int precision = 8;
   //cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&input_file, "-i", "--input",
                  "Input file to use.");

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }

   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // Create new JOTSDriver
   JOTSDriver* driver = new JOTSDriver(input_file, myid);

   // Run driver
   driver->Run();

   // Delete driver
   delete driver;

   return 0;
}