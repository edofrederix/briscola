#ifndef partialPivLUEigenSolver_H
#define partialPivLUEigenSolver_H

#include "EigenSolver.H"
#include <Eigen/Dense>

namespace Foam
{

namespace briscola
{

namespace fv
{

class partialPivLUEigenSolver
:
    public EigenSolver
{
    // Solver types
    typedef ::Eigen::PartialPivLU<EigenRhsType> EigenSolverType;

    // Solver pointer
    autoPtr<EigenSolverType> solverPtr_;

public:

    // Runtime type information
    TypeName("partialPivLU");

    // Construct
    partialPivLUEigenSolver();

    // Copy constructor
    partialPivLUEigenSolver(const partialPivLUEigenSolver&);

    // Destructor
    virtual ~partialPivLUEigenSolver();

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
    #include "partialPivLUEigenSolver.C"
#endif

#endif
