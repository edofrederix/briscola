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
    typedef ::Eigen::UmfPackLU<matrixType> solverType;

    // Solver pointer
    autoPtr<solverType> solverPtr_;

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
    virtual void prepare(const matrixType& A);

    // Solve given a right-hand side
    virtual void solve(rhsType& x, const rhsType& b) const;
};

}

}

}

#ifdef NoRepository
    #include "UmfPackEigenSolver.C"
#endif

#endif
