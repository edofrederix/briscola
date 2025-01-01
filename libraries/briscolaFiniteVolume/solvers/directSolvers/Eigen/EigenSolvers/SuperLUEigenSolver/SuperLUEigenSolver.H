#ifndef SuperLUEigenSolver_H
#define SuperLUEigenSolver_H

#include "EigenSolver.H"
#include <Eigen/SuperLUSupport>

namespace Foam
{

namespace briscola
{

namespace fv
{

class SuperLUEigenSolver
:
    public EigenSolver
{
    // Solver type
    typedef ::Eigen::SuperLU<EigenSolver::EigenMatrixType> EigenSolverType;

    // Solver pointer
    autoPtr<EigenSolverType> solverPtr_;

public:

    // Runtime type information
    TypeName("SuperLU");

    // Construct
    SuperLUEigenSolver();

    // Copy constructor
    SuperLUEigenSolver(const SuperLUEigenSolver&);

    // Destructor
    virtual ~SuperLUEigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const EigenMatrixType& A, const scalar tolerance);

    // Solve given a right-hand side
    virtual void solve(EigenRhsType& x, const EigenRhsType& b) const;

    // No first guess is needed
    virtual bool needGuess() const { return false; }
};

}

}

}

#ifdef NoRepository
    #include "SuperLUEigenSolver.C"
#endif

#endif
