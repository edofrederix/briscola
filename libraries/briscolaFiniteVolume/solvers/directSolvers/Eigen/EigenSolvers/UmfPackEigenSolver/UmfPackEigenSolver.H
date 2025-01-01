#ifndef UmfPackEigenSolver_H
#define UmfPackEigenSolver_H

#include "EigenSolver.H"
#include <Eigen/UmfPackSupport>

namespace Foam
{

namespace briscola
{

namespace fv
{

class UmfPackEigenSolver
:
    public EigenSolver
{
    // Solver type
    typedef ::Eigen::UmfPackLU<EigenMatrixType> EigenSolverType;

    // Solver pointer
    autoPtr<EigenSolverType> solverPtr_;

public:

    // Runtime type information
    TypeName("UmfPack");

    // Construct
    UmfPackEigenSolver();

    // Copy constructor
    UmfPackEigenSolver(const UmfPackEigenSolver&);

    // Destructor
    virtual ~UmfPackEigenSolver();

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
    #include "UmfPackEigenSolver.C"
#endif

#endif
