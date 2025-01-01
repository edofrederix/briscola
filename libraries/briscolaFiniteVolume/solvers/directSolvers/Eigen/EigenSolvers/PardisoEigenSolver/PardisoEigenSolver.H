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
    typedef ::Eigen::PardisoLU<EigenMatrixType> EigenSolverType;

    // Solver pointer
    autoPtr<EigenGSolverType> solverPtr_;

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
    #include "PardisoEigenSolver.C"
#endif

#endif
