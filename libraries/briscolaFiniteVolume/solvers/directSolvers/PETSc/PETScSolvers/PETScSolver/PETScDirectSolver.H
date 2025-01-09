#ifndef PETScDirectSolver_H
#define PETScDirectSolver_H

#include "PETScSolverBase.H"
#include "Switch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class PETScDirectSolver
:
    public PETScSolverBase
{
protected:

    // Factorization
    autoPtr<Mat> factor_;

public:

    // Construct from dictionary
    PETScDirectSolver(const dictionary& dict);

    // Copy constructor
    PETScDirectSolver(const PETScDirectSolver&);

    // Destructor
    virtual ~PETScDirectSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Solve given a right-hand side
    virtual void solve(Vec& x, const Vec& rhs, const MPI_Comm& comm);
};

}

}

}

#ifdef NoRepository
    #include "PETScDirectSolver.C"
#endif

#endif
