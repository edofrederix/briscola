#ifndef BiCGSTABEigenSolver_H
#define BiCGSTABEigenSolver_H

#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class BiCGSTABEigenSolver
:
    public EigenSolver
{
    // Solver types
    typedef ::Eigen::BiCGSTAB<EigenMatrixType> EigenSolverType;

    // Solver pointer
    autoPtr<EigenSolverType> solverPtr_;

public:

    // Runtime type information
    TypeName("BiCGSTAB");

    // Construct
    BiCGSTABEigenSolver();

    // Copy constructor
    BiCGSTABEigenSolver(const BiCGSTABEigenSolver&);

    // Destructor
    virtual ~BiCGSTABEigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const EigenMatrixType& A, const scalar tolerance);

    // Solve given a right-hand side
    virtual void solve(EigenRhsType& x, const EigenRhsType& b) const;

    // A first guess is needed
    virtual bool needGuess() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "BiCGSTABEigenSolver.C"
#endif

#endif
