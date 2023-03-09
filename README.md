# JOTS: MFEM Thermal Solver w/ preCICE Coupling

The MFEM unsteady thermal conduction example generalized and updated for preCICE.

## Contents
<!-- toc orderedList:0 -->

- [Dependencies](#dependencies)
- [Building JOTS](#building-jots)
        - [MFEM](#MFEM)
        - [Boost](#Boost)
        - [preCICE](#preCICE)
- [Running Simulations](#running-simulations)
        - [JOTS Configuration File](#su2-configuration-file)
        - [Running JOTS](#running-the-adapted-su2-executable)

<!-- tocstop -->

## Dependencies

### MFEM
MFEM is required - used to actually complete the FEA.

### Boost
Boost is required - used to parse input files

### preCICE
preCICE is currently required but ideally optional - used for coupling.


## Building JOTS

Dependencies include both preCICE and MFEM. The working version of this was designed using preCICE v2.5 and MFEM v4.5. When you install them, make sure you enable MPI and have the following environment variables set:


It is very important that the MPI wrappers are used to compile everything. After installing preCICE, you should have added to your bashrc the path to the preCICE install to `CMAKE_PREFIX_PATH`. It is advised that you do the same for MFEM, adding the install path to `CMAKE_PREFIX_PATH`. If you do this, and have made sure that the above environment compiler variables are set, then all that must be done is:

        mkdir build
        cd build
        cmake ..
        make

And that's it! Alternatively, you can add in arguments to preCICE and MFEM install directories if you do not have `CMAKE_PREFIX_PATH` set:

        mkdir build
        cd build
        cmake -DPRECICE_DIR=/path/to/precice/install -DMFEM_DIR=/path/to/mfem/install ..
        make

After successfully making the executable, to use the JOTS executable globally, add it to your path:

        export PATH=/path/to/jots/build:$PATH

Note that presently, JOTS requires that preCICE and MFEM headers are included and libraries linked. This may be updated in the future.

## Configuring Simulations

JOTS receives an input file .ini file. It follows a super simple format:




### Coding Plan (to be removed):

To keep as object-oriented and as conventionally-followed as possible, will set up operator class separately with .cpp and .hpp file (good practice)

Will keep JOTS as the main solver loop thing. Will likely set up an input file parser class. And maybe a restart file outputter and inputter class too. Then can easily implement additional thermal conductivity models and all that. 