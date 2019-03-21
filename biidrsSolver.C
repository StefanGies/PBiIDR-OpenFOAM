/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
Software designed to work with OpenFOAM-v1812
Software not part of official OpenFoam - Release.

Description
    Preconditioned BiIDR(S) solver with run-time selectable
    preconditioning and size of Sonneveld-Space S

Author
    Stefan Gies TU-Braunschweig

Solver is based on work by MARTIN B. VAN GIJZEN and PETER SONNEVELD, Delft University of Technology
"An Elegant IDR(s) Variant that Efficiently Exploits Biorthogonality Properties"

!BEFORE COMPILATION!
inlcude -$(CPP_DIRECTIVE) under EXE_INC in options makefile
run
> export CPP_DIRECTIVE='-std=c++11' 
or add to as compiler option in ~/foam/foam-extend-4.0/etc/bashrc
to enable c++11 compiler support
\*---------------------------------------------------------------------------*/

#include "biidrsSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(biidrsSolver, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<biidrsSolver>
        addPBiIDRAsymMatrixConstructorToTable_;
}

// conversion funtion OpenFoam Scalarfiald to Eigen C++ vector
void Foam::biidrsSolver::field2Eigen(scalarField& field, Eigen::VectorXd& out) const
 { 
     for (int l = 0; l < field.size(); l++)
     {
         out(l) = field[l];
     }
 }

// conversion funtion Eigen C++ vector to OpenFOAM Field
void Foam::biidrsSolver::Eigen2field( Eigen::VectorXd& out, scalarField& field) const
 {
     
     for (int l = 0; l < field.size(); l++)
     {
         //out(l) = field[l];
         field[l] = out(l);
     }
 }

// conversion funtion FieldField<Field, scalar> to Eigen::MatrixXd for quadratic matrices
void Foam::biidrsSolver::fieldfield2Eigen(FieldField<Field, scalar>& A, Eigen::MatrixXd& out, int xdim) const
 {

    for(int l = 0; l < xdim; l++)
    {
         for (int k = 0; k < xdim; k++)
         {
             out(l,k) = A[l][k];
         }
     }
 }



//- Construct from matrix and solver data stream
Foam::biidrsSolver::biidrsSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    sDims_(readLabel(solverControls.lookup("sDimensions"))),
    kappa_(readScalar(solverControls.lookup("angle"))),
    resprint_(readLabel(solverControls.lookup("resprint"))),
    subspace_(word(solverControls.lookup("subSpace")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * /
// manual loop unrolling? does it lead to performance increase? probs not

Foam::solverPerformance Foam::biidrsSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const

{
      // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    autoPtr<lduMatrix::preconditioner> preconPtr =
    lduMatrix::preconditioner::New
    (
        *this,
        controlDict_
    );


    const label nCells = x.size();
    const scalar* __restrict__ bPtr = b.begin();

    scalar* __restrict__ xPtr = x.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalarField tempField(nCells);

    //  wA = Ax;
    matrix_.Amul(wA, x, interfaceBouCoeffs_, interfaces_, cmpt);

    // calculate normfactor
    const scalar normFactor = this->normFactor(x, b, wA, tempField);

    if (lduMatrix::debug >= 2)
    {
        Info<< "Normalisation factor = " << normFactor << endl;
    }

    // r = b - A*x_o;
    scalarField rA(nCells);
    scalar* __restrict__ rAPtr = rA.begin();

    for (register label cell = 0; cell < nCells; cell++)
    {
        rAPtr[cell] = bPtr[cell] - wAPtr[cell];
    }
    


    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (resprint_ == 1)
    {
        Info << solverPerf.initialResidual() << endl;
    }
   

   // check convergence
    if 
    (   minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
         
    // set N ---- size of system
    int S = sDims_;
    
    // set matrix P -> basis of generated Krylov Subspace
    FieldField <Field, scalar> P(S);
    forAll (P, i)
    {
        P.set(i, new scalarField(nCells, 0.0)); // inititlayse as 0.0
    }
    
    label i,j;


    Random randomA(1234567);

    // choose composition of subspace P = [r0,rand] or P = [rand]
    if(subspace_ == "r0")
    {
        // P[0] = r0 to approximate BiCGStab convergence behaviour
        scalar* __restrict__ PnullPtr = P[0].begin();

        for (register label cell = 0; cell < nCells; cell++)
        {
            PnullPtr[cell] = rAPtr[cell];
        }


        if (S > 1)
        {
            for (i = 1; i < S; i++)
            {
                    
                scalar* __restrict__ PiPtr = P[i].begin();

                // full random P[i] i > 0  
                for (register label cell = 0; cell < nCells ; cell++)
                {
                    PiPtr[cell] = randomA.GaussNormal<scalar>();
                }    
            }   
        }
    }
    else if (subspace_ == "rand") 
    {
        for (i = 0; i < S; i++)
        {
                
            scalar* __restrict__ PiPtr = P[i].begin();

            // full random P[i] i >= 0  
            for (register label cell = 0; cell < nCells ; cell++)
            {
                PiPtr[cell] = randomA.GaussNormal<scalar>();
            }    
        }   
    } 
    else 
    {
        FatalErrorInFunction
            << "PBiIDR: please select 'r0' or 'rand' for the parameter subSpace"
            << exit(FatalError);
    }     
    
    //  orthogonalisation of matrix P --> initial Shadow Space Space S in K_0 
    // might want to replace mod. Gram-Schmidt Procedure with sth more accurate
    label lk ;
    long double alpha = 0.0;

    for (lk = 1; lk < S; ++lk)
    {

        scalar* __restrict__ PlkPtr = P[lk].begin();
        
        for(j = 0; j < lk; ++j)
        {

            scalar* __restrict__ PjPtr = P[j].begin();
            double vu = gSumProd(P[lk],P[j], matrix().mesh().comm());
            double uu = gSumSqr(P[j],  matrix().mesh().comm());   

            if (!solverPerf.checkSingularity(mag(uu)))
            {
                alpha = vu / uu;
            }
            else
            {
                break;
            }

            for(register label cell = 0; cell < nCells; cell++)
            {
                PlkPtr[cell] -= alpha*PjPtr[cell];
            }

        }

      }
    
      

    //  Mu - Matrix  - Projection of residuals onto Shadow Space
    FieldField <Field, scalar> Mu(S);
    forAll (Mu, i)
    {
        Mu.set(i, new scalarField(S, 0.0)); // initialyse as 0.0
    } 

    // initialyse Mu as eye(S)
    for (i = 0; i < S; i++)
    {
        Mu[i][i] = 1.0;
    }
    
    
 
    //  omega - scaling of residual projection
    scalar omega = 1.0;


    // initialyse phi
    scalarField phi(S,0.0);

    // initialyse t = Av
    scalarField t(nCells);   
    scalar* __restrict__ tPtr = t.begin();

    //search direction v
    scalarField v(nCells);
    scalar* __restrict__ vPtr = v.begin();   
   
    //search direction pv
    scalarField pv(nCells);
    scalar* __restrict__ pvPtr = pv.begin();         
    
    //norms and scalar products initialisation
    scalar normr = 0.0;
    scalar normt = 0.0;
    scalar tr = 0.0;
    scalar rho = 0.0;
    

    // set U, containing s vectors of N in C = 0
    FieldField <Field, scalar> U(S);
    forAll (U, i)
    {
        U.set(i, new scalarField(nCells, 0.0));
    }
    
    // set G(N, s) = U(N, s) = 0;
    FieldField <Field, scalar> G = U;
    scalar beta = 1.0;

    do
    {
        label k, l;

        // phi = P * rA
        for (l = 0; l<S; l++)
        {

            phi[l] = gSumProd(P[l],rA,matrix().mesh().comm());
        }


        // iterate through s-dimensional subspace
        for (k = 0; k < S; ++k)
        {

            // check for convergence every time a new intermediate res 
            // was generated
            if (!solverPerf.checkConvergence(tolerance_, relTol_)) 
            {
                
                // pointers to elements in G_k and U_k
                scalar* __restrict__ UkPtr = U[k].begin();
                scalar* __restrict__ GkPtr = G[k].begin();
                    
                // calculate dimension of small system
                int sysdim = S-k; 
                label i, j;

                // set new small system M_s * c = phi_s to solve
                FieldField <Field, scalar> Ms(sysdim);
                forAll (Ms, i)
                {
                    Ms.set(i, new scalarField(sysdim, 0.0)); // inititlayse as 0.0
                }

                // update M and phi for solution of c = M^(-1)Phi
                for (i = k; i < S ; i++)
                {
                    for(j = k; j < S; j++)
                    {
                        Ms[i-k][j-k] = Mu[i][j];
                    }
                }
                    
                scalarField phis(sysdim);
                for (i = k ; i < S; i++)
                {
                    phis[i-k] = phi[i];
                }

                // solution vector c 
                scalarField c(sysdim,0.0);
   
                // c = Mu⁻¹*phi ; generate Eigen::matrix and vector containers
                Eigen::MatrixXd ME(sysdim, sysdim);
                Eigen::VectorXd phisE(sysdim);
                Eigen::VectorXd cE(sysdim);

                // convert Foam::Fields to Eigen
                Foam::biidrsSolver::field2Eigen(phis, phisE);
                Foam::biidrsSolver::fieldfield2Eigen(Ms, ME, sysdim);

                // use Eigen::SVD to solve small system c = Mu⁻¹*phi 
                Eigen::BDCSVD<Eigen::MatrixXd> 
                    svd(ME,  Eigen::ComputeFullU | Eigen::ComputeFullV);
                cE = svd.solve(phisE);

                // convert solution back to scalarField
                Foam::biidrsSolver::Eigen2field(cE, c);
        
                    
                scalarField ukc(nCells,0.0);
                scalar* __restrict__ ukcPtr = ukc.begin();
                scalarField gkc(nCells,0.0);
                scalar* __restrict__ gkcPtr = gkc.begin();

                // compute G[k]*c and U[k]*c;
                for (i = k; i < S; i++)
                {
                    label ck = i-k;
                    scalar* __restrict__ GiPtr = G[i].begin();
                    scalar* __restrict__ UiPtr = U[i].begin();

                    
                    for (register label cell = 0; cell < nCells; cell++)
                    {
                        gkcPtr[cell] += GiPtr[cell]*c[ck];
                        ukcPtr[cell] += UiPtr[cell]*c[ck];
                    }
                        
                }

                    
                //v = rA - G[k]*c
                for (register label cell = 0; cell < nCells; cell++)
                {
                    vPtr[cell] = rAPtr[cell] - gkcPtr[cell];
                }            

                //preconditioning operation on search direction v
                preconPtr->precondition(pv, v, cmpt);
                    
                //U[k] =  omega*pv + U[k}*c;
                for (register label cell = 0; cell < nCells; cell++)
                {
                    UkPtr[cell] = omega*pvPtr[cell] + ukcPtr[cell];
                }

                //G[K] = A*U[k]
                matrix_.Amul(G[k], U[k], interfaceBouCoeffs_, interfaces_, cmpt);


                // Gram Schmidt - Type Orthogonalisation 
                for (i = 0; i < k; i++)
                {                        
                        
                    alpha = 
                        gSumProd(P[i],G[k], matrix().mesh().comm())
                            /Mu[i][i];

                    scalar* __restrict__ GiPtr = G[i].begin();
                    scalar* __restrict__ UiPtr = U[i].begin();

                    
                    for (register label cell = 0; cell < nCells; cell++)
                    {
                        GkPtr[cell] -= alpha*GiPtr[cell];
                        UkPtr[cell] -= alpha*UiPtr[cell];
                    }
                }

                    
                // Mu = P_i*G_k
                for (i = k; i < S; i++)
                {
                    Mu[i][k] = gSumProd(P[i],G[k]);
                }
                    
                // beta = phi[k]/Mu[k][k]    
                if(Mu[k][k] != 0.0)
                {
                    beta = phi[k]/Mu[k][k];

                }
                else
                {
                    Info << "Breakdown of BiIDR(" << S 
                        <<") at stage "<< k << endl; 
                    solverPerf.finalResidual() =
                            gSumMag(rA, matrix().mesh().comm())/normFactor;  
                    solverPerf.nIterations()++;
                    return solverPerf;   
                }

                // r = r - beta*G_k
                // x = x + beta*U_k
                for (register label cell = 0; cell < nCells; cell++)
                {
                    rAPtr[cell] -= beta*GkPtr[cell];
                    xPtr[cell] += beta*UkPtr[cell];
                }

                    
            }
            else
            {
                // if norm(r) < TOL --> end solver with intermediate res as final
                solverPerf.finalResidual() 
                    = gSumMag(rA, matrix().mesh().comm())/normFactor;
                solverPerf.nIterations()++;
                return solverPerf;
            }     

            label i; 
            // update RHS of small system                
            if (k < S - 1)
            {   
                for (i = 0; i < S ; i++)
                {  
                    if (i > k)
                    {
                        phi[i] -= beta*Mu[i][k];
                    }
                    else
                    {
                        phi[i] = 0.0;
                    }
                }
            }

        // finish loop over intersection Gi and Sonneveld Space           
        }
                
        //preconditioning operation
        preconPtr->precondition(v, rA, cmpt);
        // t = Av
        matrix_.Amul(t, v, interfaceBouCoeffs_, interfaces_, cmpt);

        // calculate L1 norms of rA and t
        normr = gSumMag(rA, matrix().mesh().comm());
        normt = gSumMag(t, matrix().mesh().comm());

        // choose omega to minimise ||r_{i+1}|| = ||r_{i} - omega{i+1}*t||
        omega = normr/normt;

        
        // execute omega stablilisation with angle = kappa
        // see Sonneveld et al for Details
        if (kappa_ != 0.0)
        {
            scalar l2normr = Foam::sqrt(gSumSqr(rA, matrix().mesh().comm()));
            scalar l2normt = Foam::sqrt(gSumSqr(t, matrix().mesh().comm()));

            tr = gSumProd(t,rA, matrix().mesh().comm());
            scalar normtr = abs(tr);
                
            rho = (normtr)/(l2normr*l2normt);
                
            if (solverPerf.checkSingularity(mag(rho)))
            {
                break;
            }

            if (rho < kappa_)
            {
                omega = omega*(kappa_/rho); 
            }
                
        }

           
        // projection onto lower dim Krylov Space
        // (I-omegaA)(G(n+1) intersection S) --> update residual and solution
        for (register label cell = 0; cell < nCells; cell++)
        {
            rAPtr[cell] -= omega*tPtr[cell];
            xPtr[cell] += omega*vPtr[cell]; 
        }    
                      
    
        // evaluate final residual and increase iteration count
        solverPerf.finalResidual() = gSumMag(rA)/normFactor;
        //solverPerf.finalResidual() = sqrt(gSumProd(rA,rA))/normFactor;
        solverPerf.nIterations()++;

        if(resprint_ == 1)
        {
            Info << solverPerf.finalResidual() << endl;
        }

    } while (   
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
         );
    }

    // Recommend adjustments if PBiIDR fails to converge
    if (solverPerf.nIterations() > max(defaultMaxIter_, maxIter_))
    {
        FatalErrorInFunction
            << "PBiIDR has failed to converge within the maximum number"
               " of iterations " << max(defaultMaxIter_, maxIter_) << nl
            << "    Please try reducing the subspace-dimensionality sDims or adjusting the 'angle' parameter."
            << exit(FatalError);
    }

    matrix().setResidualField(rA, fieldName_, false);

    return solverPerf;
}


// ************************************************************************* //
