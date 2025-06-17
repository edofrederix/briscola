#ifndef PETScExternalSolver_H
#define PETScExternalSolver_H

#include "PETScSolverBase.H"
#include "Switch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class PETScExternalSolver
:
    public PETScSolverBase
{
protected:

    // Factorization
    autoPtr<Mat> factor_;

public:

    // Construct from dictionary
    PETScExternalSolver(const dictionary& dict);

    // Copy constructor
    PETScExternalSolver(const PETScExternalSolver&);

    // Destructor
    virtual ~PETScExternalSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Solve given a right-hand side
    virtual void solve(Vec& x, const Vec& rhs, const MPI_Comm& comm);
};

}

}

}

#ifdef NoRepository
    #include "PETScExternalSolver.C"
#endif

#endif
