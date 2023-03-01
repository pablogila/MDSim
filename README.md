# Molecular Dynamics Simulations

The program MDS_init.f90 generates the initial positions for a 3D cube of atoms, and MDS_EOM.f90 computes the Equations Of Motion of said cube via Verlet algorithm, taking Lennard-Jones interactions into account.

Download the last version from:\
https://github.com/pablogila/Molecular_Dynamics_Sim

Since outputs can become too large, it is recommended to clean everything with RepoCleaner after saving the relevant data:\
https://github.com/pablogila/RepoCleaner


## Instructions

Compile the program and execute as follows:

First, we assign the initial positions, creataing the initial files for the cube\
`gfortran tipos.f90 MDS_init.f90`\
`./a.out`

Next we execute the Equations Of Motion\
`gfortran tipos.f90 MDS_EOM.f90`\
`./a.out`


## Visualization

Visualize initial cube\
`xmakemol -f fort.3`

Plot histogram of velocities\
`gnuplot`\
`> plot 'fort.7' u 1:2 w l`

Plot histogram of positions\
`gnuplot`\
`> plot 'fort.8' u 1:2 w l`

Plot the kinetic, potential and total energies\
`gnuplot`\
`> plot "fort.10" u 2:3 w l, "fort.10" u 2:4 w l, "fort.10" u 2:5 w l`

Visualize the simulation\
`xmakemol -f fort.11`


## Files created

fort.3  - initial positions in xyz, for xmakemol\
fort.4  - initial positions\
fort.7  - histogram of velocities\
fort.8  - histogram of positions\
fort.9  - initial velocities\
fort.10 - energies at each time step\
fort.11 - positions at each time step
