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
    typedef ::Eigen::SuperLU<matrixType> solverType;

    // Solver pointer
    autoPtr<solverType> solverPtr_;

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
    virtual void prepare(const matrixType& A);

    // Solve given a right-hand side
    virtual void solve(rhsType& x, const rhsType& b) const;
};

}

}

}

#ifdef NoRepository
    #include "SuperLUEigenSolver.C"
#endif

#endif
