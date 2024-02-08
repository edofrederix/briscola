#ifndef PardisoEigenSolver_H
#define PardisoEigenSolver_H

#include "EigenSolver.H"
#include <Eigen/PardisoSupport>
#include <Eigen/Core>

namespace Foam
{

namespace briscola
{

namespace fv
{

class PardisoEigenSolver
:
    public EigenSolver
{
    // Solver types
    typedef ::Eigen::PardisoLU<matrixType> gSolverType;
    typedef ::Eigen::PardisoLDLT<matrixType> sSolverType;

    // Solver pointers
    autoPtr<gSolverType> gSolverPtr_;
    autoPtr<sSolverType> sSolverPtr_;

public:

    // Runtime type information
    TypeName("Pardiso");

    // Construct
    PardisoEigenSolver();

    // Copy constructor
    PardisoEigenSolver(const PardisoEigenSolver&);

    // Destructor
    virtual ~PardisoEigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const matrixType& A);

    // Solve given a right-hand side
    virtual void solve(rhsType& x, const rhsType& b) const;
};

}

}

}

#ifdef NoRepository
    #include "PardisoEigenSolver.C"
#endif

#endif
