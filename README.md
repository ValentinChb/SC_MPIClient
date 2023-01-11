# SC_MPIClient

This is a MPI bridge enabling alternative supercontroller implementations for NREL's FAST.farm (https://github.com/OpenFAST/openfast). 
It establishes a two-way real-time connection between on the one side a controller run in an external parallel environment (e.g. The Mathworks' Matlab/Simulink), and on the other side OpenFAST instances for each turbine in FAST.Farm, through the bladed-style controller dll loaded in ServoDyn.

At each turbine-level (ServoDyn) timestep, the bladed-style controller array (avrSWAP) is exchanged between the client (i.e. the turbine side, provided in this repository) and the server (i.e. the farm side).
Parallel execution of OpenFAST instances within one farm-level (FAST.Farm) timestep is supported (using the OMP version of FAST.Farm). 

It may be linked to an existing bladed-style turbine controller implementation (DTUWEC or ROSCO) to provide turbine-level controls (generator torque and blade pitch angles), while the external controller provides farm-level controls (yaw angles and power setpoints). 
If not, both farm- and turbine-level controls are read from the external controller. 

# Installation
- Install a MPI distribution (Microsoft MPI has been used in this project) and run Build_lmsmpi.bat (adapt if necessary) to generate linkable MPI dependencies
- Use CMAKE to generate Makefile, with *src* as source folder and *custom-build* as build folder. Specific CMAKE BOOL options are:
  - CMAKE_LINK2DTUWEC: links to a modified version of the DTU Wind Energy Controller (https://github.com/ValentinChb/DTUWEC) for turbine-level controls, implementing an active derating functionality
  - CMAKE_LINK2ROSCO: links to the ROSCO controller (https://github.com/NREL/ROSCO) as an alternative/complement for turbine-level controls
  - CMAKE_LINK2MPI: if deactivated, no MPI bridge with external controller will be made. This may be used for single-turbine simulations to avoid unnecessary dependencies, or for debugging
- Make sure you have the following dependencies (update if actual) in custom-build/src: mpi.mod and libmsmpi.a (if using other mpi distribution, update name in CMAKE), dtu_we_controller_bladed.mod and libDTUWEC4SC.a (if using DTUWEC), rosco.mod, rosco_types.mod and libROSCO.a (if using ROSCO). Note that linking must be static, as dynamic linking (loading DTUWEC or ROSCO dll at runtime as dependency) is challenging to make thread-safe.
- Make project. The library will have a different name depending on the CMAKE options above.


# Use
- FAST.Farm input files must satisfy the following architecture: <Root folder containing .fstf file> / OpenFAST / T< iturb > / <.fst file for turbine iturb>
- The input file SC_input.dat must be placed in the OpenFAST folder, with the following inputs: 
  - UseSC (0 or 1) use MPI communication or not. 0 has similar effect as CMAKE_LINK2MPI=FALSE at build stage, except that the library is linked to depednencies though not using them
  - Number of turbines
  - Farm-level timestep
  - Path of the MPI communication shared file, used to sync communication information between the MPI client and server (you should not edit and should not need to open this file, just make sure its parent folder exists)
  - Mod_AmbWind in FAST.Farm (for advanced use, not supported here)
  - Nseeds for running multiple realizations of turbulent wind field (for advanced use, not supported here)
