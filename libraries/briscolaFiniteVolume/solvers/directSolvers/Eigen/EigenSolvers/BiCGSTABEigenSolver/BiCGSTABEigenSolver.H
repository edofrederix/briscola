#ifndef BiCGSTABEigenSolver_H
#define BiCGSTABEigenSolver_H

#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef ::Eigen::BiCGSTAB<EigenSolverBase::RowMajorMatrixType>
    BiCGSTABSolver;

class BiCGSTABEigenSolver
:
    public EigenSolver<BiCGSTABSolver,::Eigen::RowMajor>
{
    // Maximum number of iterations (defaults to 1000)
    label maxIter_;

    // Tolerance (defaults -1, meaning that the default one is used)
    scalar tolerance_;

public:

    // Runtime type information
    TypeName("BiCGSTAB");

    // Construct from dictionary
    BiCGSTABEigenSolver(const dictionary& dict);

    // Copy constructor
    BiCGSTABEigenSolver(const BiCGSTABEigenSolver&);

    // Destructor
    virtual ~BiCGSTABEigenSolver();

    // Solve given a right-hand side
    virtual void solve(RhsType& x, const RhsType& b);

    // First guess is needed
    virtual bool needGuess() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "BiCGSTABEigenSolver.C"
#endif

#endif
