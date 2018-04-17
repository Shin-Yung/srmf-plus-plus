About
-------
A C++ program for diagonalizing quadratic tight binding Hamiltonian.


Version
-------
1.0.0 - June 2017.

Build
------------
To build, first make a copy of the file 'make_options.template' in the 
project root directory and rename it as 'make_options.mk'. Edit this file 
to set your variables. There are two dependencies - Boost and Eigen C++. 
The C++ compilers need to support C++11 features.  

Usage
-----
To run the application (say, ./a.out), do 

	./a.out [OPTIONS] [FILE]  

The program accepts a few optional arguments OPTIONS. Type 

	./a.out -help

to see the list of available options. The FILE argument, if present, 
is taken as the input filename, else a default filename 'input.parm' is 
looked for in the current directory.  
