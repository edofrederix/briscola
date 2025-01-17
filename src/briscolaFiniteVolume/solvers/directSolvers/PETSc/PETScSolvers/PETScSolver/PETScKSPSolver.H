#ifndef PETScKSPSolver_H
#define PETScKSPSolver_H

#include "PETScSolverBase.H"
#include "Switch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class PETScKSPSolver
:
    public PETScSolverBase
{
protected:

    // Solver pointer
    autoPtr<KSP> solverPtr_;

    // Absolute tolerance (defaults to zero)
    scalar tolerance_;

    // Relative tolerance (defaults to 1e-5)
    scalar relTol_;

    // Maximum number of iterations (defaults to 1000)
    label maxIter_;

    // Print statistics (defaults to false)
    Switch printStats_;

public:

    // Construct from dictionary
    PETScKSPSolver(const dictionary& dict);

    // Copy constructor
    PETScKSPSolver(const PETScKSPSolver&);

    // Destructor
    virtual ~PETScKSPSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Solve given a right-hand side
    virtual void solve(Vec& x, const Vec& rhs, const MPI_Comm& comm);
};

}

}

}

#ifdef NoRepository
    #include "PETScKSPSolver.C"
#endif

#endif
