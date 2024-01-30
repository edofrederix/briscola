#ifndef EigenSolver_H
#define EigenSolver_H

#include "dictionary.H"
#include "autoPtr.H"

#include <Eigen/Sparse>

namespace Foam
{

namespace briscola
{

namespace fv
{

class EigenSolver
{
public:

    // Eigen types
    typedef ::Eigen::SparseMatrix<double> matrixType;
    typedef ::Eigen::MatrixXd rhsType;

protected:

    // Check if the given sparse matrix is symmetric
    bool isSymmetric(const matrixType&) const;

    // Check if the given sparse matrix is positive definite
    bool isPositiveDefinite(const matrixType&) const;

public:

    // Runtime type information
    TypeName("EigenSolver");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        EigenSolver,
        dictionary,
        (),
        ()
    );

    // Construct
    EigenSolver();

    EigenSolver(const EigenSolver&);

    // Selector from solver type
    static autoPtr<EigenSolver> New(const word);

    // Destructor
    virtual ~EigenSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const matrixType& A) = 0;

    // Solve given a right-hand side
    virtual void solve(rhsType& x, const rhsType& b) const = 0;
};

}

}

}

#endif
