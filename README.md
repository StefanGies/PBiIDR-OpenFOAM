# PBiIDR-OpenFOAM
Preconditioned biorthogonal IDR(s) Solver for OpenFoam 

Solver uses Eigen3 C++ Library for Dense Matrix Algebra
  current version uses Block-Divide and Conquer SVD
  
compiled with Options: -O3 -march=native -march=native -mavx2 (if available) // change -march and -mtune accordingly
DO NOT compile with -Ofast or -ffast-math > any option that changes floating point behaviour
to a lower precision
-frounding-maths is possible, but might have detrimental effects on performance in not yet tested cases


Solver is an implementation of work by MARTIN B. VAN GIJZEN and PETER SONNEVELD, Delft University of Technology
"An Elegant IDR(s) Variant that Efficiently Exploits Biorthogonality Properties"