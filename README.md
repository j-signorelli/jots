# JOTS v2.0: MFEM-Based Thermal Solver w/ preCICE Coupling
JOTS is a heat transfer solver designed for either standalone heat transfer numerical simulations or conjugate heat transfer analyses coupled using preCICE.

## Contents
<!-- toc orderedList:0 -->

- [Dependencies](#dependencies)
    - [MFEM](#MFEM)
    - [Boost](#Boost)
    - [preCICE](#preCICE)
- [Building JOTS](#building-jots)
- [Running Simulations](#running-simulations)

<!-- tocstop -->

## Dependencies

### MFEM
JOTS requires MFEM v4.7. Detailed build instructions for MFEM can be found here: https://mfem.org/building/.

Note that you *must* build MFEM with `MFEM_USE_MPI=YES`. If you want compressed VisIt/Restart files as opposed to ASCII, build MFEM with `MFEM_USE_ZLIB=YES`.

### Boost
Boost headers are used to parse input files. JOTS was developed using v1.82.0.

### preCICE
preCICE is used for coupling JOTS with fluid solvers. JOTS was developed using v2.5. Detailed instructions for building preCICE can be found here: https://precice.org/installation-source-preparation.html.

## Building JOTS

After building + installing the above dependencies, you may need to add the install locations for each of them to `CMAKE_PREFIX_PATH`. Then:

        export CXX=$(which mpicxx)
        mkdir build
        cd build
        cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ ..
        make -j
        make install

Alternatively, you can add in arguments to install directories if you do not have `CMAKE_PREFIX_PATH` set with them, like:

        cmake -DPRECICE_DIR=/path/to/precice/install -DMFEM_DIR=/path/to/mfem/install ..

After successfully making the executable, to use the JOTS executable globally, add to bashrc:

        export PATH=/path/to/jots/install/bin/:$PATH

Note that the pw_mfem_mesh_fix.py script is also copied and pasted in the bin folder of the install directory.

To enable regression testing using ctest, include `-DENABLE_TESTING=ON` in the `cmake` command.

## Running Simulations

JOTS using a .ini file as an input file. All boundary conditions and settings can be found in the template config files.

To run JOTS with 10 procs:

        mpirun -np 10 jots -i config_file.ini