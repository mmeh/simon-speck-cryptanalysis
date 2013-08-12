SIMON/SPECK cryptanalysis code
==============================

Cryptanalysis code for the SIMON and SPECK families of block ciphers

Authors
-------
This SIMON implementation and all cryptanalytic routines found in "simon.h" and "simon.cpp" were written by 

* Martin M. Lauridsen   (mmeh@dtu.dk) and 
* Hoda A. Alkhzaimi		(hoalk@dtu.dk)

DTU Compute  
Section for Cryptology  
Department of Applied Mathematics and Computer Science  
Matematiktorvet, building 324

Compilation
-----------
To compile the code, use the attached Makefile. You will see, that compilation flags include the C++11 standard and the OpenMP library ("omp.h"), so make sure you are using a version of gcc that supports C++11, and that you have the OpenMP library installed (or alternatively remove the parts of the code that uses OpenMP).

Documentation
-------------
We have not been very good at documenting the code, but in "simon.h", the different methods are categorized by their role in the different cryptanalytic approaches (differential, impossible differential, etc.).
