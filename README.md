# JOTS: MFEM Thermal Solver w/ preCICE Coupling

The MFEM unsteady thermal conduction example generalized and updated for preCICE.

## Contents
<!-- toc orderedList:0 -->

- [Building JOTS](#building-jots)
- [Configuring Simulations](#configuring-simulations)
    - [SU2](#su2)
    - [preCICE](#precice)
    - [Adapter](#adapter)
- [Running Simulations](#running-simulations)
    - [SU2 Configuration File](#su2-configuration-file)
    - [Running the Adapted SU2 Executable](#running-the-adapted-su2-executable)
- [References](#references)

<!-- tocstop -->

## Building JOTS

Dependencies include both preCICE and MFEM. The working version of this was designed using preCICE v2.5 and MFEM v4.5. When you install them, make sure you enable MPI and have the following environment variables set:

        export CC=$(which mpicc)
        export CXX=$(which mpicxx)
        export FC=$(which mpif90)

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

JOTS ca