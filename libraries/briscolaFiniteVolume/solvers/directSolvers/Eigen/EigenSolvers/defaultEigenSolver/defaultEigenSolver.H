#ifndef defaultEigenSolver_H
#define defaultEigenSolver_H

#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class defaultEigenSolver
:
    public EigenSolver
{
    // Solver types
    typedef ::Eigen::SparseLU<matrixType> gSolverType;
    typedef ::Eigen::SimplicialLDLT<matrixType> sSolverType;

    // Solver pointers
    autoPtr<gSolverType> gSolverPtr_;
    autoPtr<sSolverType> sSolverPtr_;

public:

    // Runtime type information
    TypeName("default");

    // Construct
    defaultEigenSolver();

    // Copy constructor
    defaultEigenSolver(const defaultEigenSolver&);

    // Destructor
    virtual ~defaultEigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const matrixType& A);

    // Solve given a right-hand side
    virtual void solve(rhsType& x, const rhsType& b) const;
};

}

}

}

#ifdef NoRepository
    #include "defaultEigenSolver.C"
#endif

#endif
