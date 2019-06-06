# SSS-VSSM sample code

https://github.com/hessfran/SSS-VSSM

Author: F. Hess (dr.franziska.hess@gmail.com) www.hessgroup-catalysis.de


These programs are used to showcase and test the SSS (Supercluster, Subtraction, Supersite) algorithms. These algorithms are intended as a modification to the Variable Step Size Method (VSSM) to enable more efficient surface Kinetic Monte Carlo (KMC) simulations with Cluster Expansion (CE)-based lateral interactions, but there may be other possible applications.

The CE model assumes a rectangular lattice (C2v symmetry) such as the RuO2(110) surface with four species + vacant site (five in total).

More detailed instructions and information are found in each of the program files.



Feel free to contact me if there is any problem and please notify me of any bugs you find.


Compilation
-----------

The programs were tested for ifort 11.1 and ifort 18.0.1  (no flags required).

Only supercluster.f90 can be used with gfortran and requires the -ffree-line-length-200 flag:

gfortran -o supercluster -ffree-line-length-200 supercluster.f90