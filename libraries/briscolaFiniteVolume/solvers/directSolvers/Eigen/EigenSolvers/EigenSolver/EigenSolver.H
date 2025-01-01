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

    // Eigen matrix type with row-major ordering
    typedef typename ::Eigen::SparseMatrix<double> EigenMatrixType;

    // Eigen right-hand side type
    typedef typename ::Eigen::MatrixXd EigenRhsType;

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
    virtual void prepare(const EigenMatrixType& A, const scalar tolerance) = 0;

    // Solve given a right-hand side
    virtual void solve(EigenRhsType& x, const EigenRhsType& b) const = 0;

    // Specifiy if a first guess is needed
    virtual bool needGuess() const = 0;
};

}

}

}

#endif
