# Molecular Dynamics Simulations
This program calculates a 3D cube via Lennard-Jones potential, and computes the Equations Of Motion via Verlet algorithm.

# Instructions:
Compile the program and execute as follows:
* `gfortran tipos.f90 init-positions.f90`
* `./a.out`

Now that we have created the initial files for the cube, we execute the Equations Of Motion:\\
* `gfortran tipos.f90 EOM.f90`
* `./a.out`

And visualize the simulation:
* `xmakemol -f fort.11`


To plot kinetic energies:
* `gnuplot`
* `> plot "fort.10" u 2:3 w l, "fort.10" u 2:4 w l, "fort.10" u 2:5 w l`
