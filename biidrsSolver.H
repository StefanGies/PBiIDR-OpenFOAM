/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This software is inteded for OpenFOAM.

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This Software is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    biidrsSolver

Description
    Preconditioned Bi-Orthogonal Inducded Dimension Reduction s-stage solver with run-time 	selectable preconditioning.

Author
    Stefan Gies, TU-Braunschweig

SourceFiles
    biidrsSolver.C

Solver is based on work by MARTIN B. VAN GIJZEN and PETER SONNEVELD
, Delft University of Technology

\*---------------------------------------------------------------------------*/

#ifndef biidrsSolver_H
#define biidrsSolver_H

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include "Random.H"
#include <omp.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>



using namespace std;

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class biIDR(s)Solver Declaration
\*---------------------------------------------------------------------------*/

class biidrsSolver
:
    public lduMatrix::solver
{
    // Private data

        label sDims_;
        scalar kappa_;
        label resprint_; 
        word subspace_;


    // Private Member Functions

    //- Disallow default bitwise copy construct
    biidrsSolver(const biidrsSolver&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const biidrsSolver&) = delete;

    void field2Eigen(
        Field<scalar>& field, Eigen::VectorXd& out) const;

    void Eigen2field(
        Eigen::VectorXd& out, scalarField& field) const;

    void fieldfield2Eigen(
        FieldField<Field, scalar>& A, Eigen::MatrixXd& out, int xdim) const;



public:

    //- Runtime type information
    TypeName("PBiIDR");


    // Constructors
    //- Construct from matrix components and solver data stream
      biidrsSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    // Destructor

    virtual ~biidrsSolver()
       {}


    // Member Function

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
