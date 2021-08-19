# PBiIDR-OpenFOAM
Preconditioned biorthogonal IDR(s) Solver for OpenFoam (tested for versions v1812, v1912 and v2012)i

Solver uses Eigen3 C++ Library for Dense Matrix Algebra
  current version uses Block-Divide and Conquer SVD
  
compiled with Options: -O3 -march=native -march=native -mavx2 (if available) // change -march and -mtune accordingly
DO NOT compile with -Ofast or -ffast-math > any option that changes floating point behaviour
to a lower precision
-frounding-maths is possible, but might have detrimental effects on performance in not yet tested cases


Solver is an implementation of work by MARTIN B. VAN GIJZEN and PETER SONNEVELD, Delft University of Technology
"An Elegant IDR(s) Variant that Efficiently Exploits Biorthogonality Properties"

USAGE GUIDE:
As with any other linear algebra solver in OpenFOAM the solver is called from the fvSolution dictionary
Example:
```
U
  {
      solver          PBiIDR;
      preconditioner  none;
      tolerance       1e-05;
      relTol          0;
      sDimensions	    2;
      subSpace	      rand;   
      angle		        0.7;  
        resprint	      0;	  
  }
```
  sDimensions:  dimensionality of Sonneveld subspace
  subSpace:     first vector to build Sonneveld space with > span{r0/rand, rand, ..., rand}, choose 'r0' or 'rand'; if sDimension > 1 -> coloum vectors (full) random numbers (Gaussian) distribution
  angle:        adjustment of omega based on angle between prior residual and search direction; if angle == 0.0, omega is adjusted to minimise residual norm
  resprint:     print residuals  at every internal iteration in std output (debugging only)


Disclaimer:
This software is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
