#ifndef EigenSolver_H
#define EigenSolver_H

#include "EigenSolverBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SolverType, int Order>
class EigenSolver
:
    public EigenSolverBase
{
public:

    // Eigen matrix type
    typedef typename ::Eigen::SparseMatrix<double, Order> MatrixType;

protected:

    // Solver pointer
    autoPtr<SolverType> solverPtr_;

public:

    // Construct from dictionary
    EigenSolver(const dictionary& dict);

    // Copy constructor
    EigenSolver(const EigenSolver<SolverType,Order>&);

    // Destructor
    virtual ~EigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const MatrixType& A);

    // Solve given a right-hand side
    virtual void solve(RhsType& x, const RhsType& b);

    // Storage order
    virtual int order() const { return Order; }
};

}

}

}

#ifdef NoRepository
    #include "EigenSolver.C"
#endif

#endif
