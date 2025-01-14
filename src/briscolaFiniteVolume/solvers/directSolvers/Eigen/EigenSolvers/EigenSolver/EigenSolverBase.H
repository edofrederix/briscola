#ifndef EigenSolverBase_H
#define EigenSolverBase_H

#include "dictionary.H"
#include "autoPtr.H"
#include <Eigen/Sparse>

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class EigenLinearSystem;

class EigenSolverBase
{
protected:

    // Copy of the dictionary
    dictionary dict_;

public:

    // Eigen types

    typedef typename ::Eigen::SparseMatrix<double, ::Eigen::RowMajor>
        RowMajorMatrixType;
    typedef typename ::Eigen::SparseMatrix<double, ::Eigen::ColMajor>
        ColMajorMatrixType;

    typedef typename ::Eigen::MatrixXd RhsType;

    // Runtime type information
    TypeName("EigenSolverBase");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        EigenSolverBase,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Construct from dictionary
    EigenSolverBase(const dictionary& dict);

    // Copy constructor
    EigenSolverBase(const EigenSolverBase&);

    // Selector from solver type and coefficient dictionary
    static autoPtr<EigenSolverBase> New(const word, const dictionary& dict);

    // Destructor
    virtual ~EigenSolverBase();

    // Prepare from a sparse matrix
    virtual void prepare(const RowMajorMatrixType& A);
    virtual void prepare(const ColMajorMatrixType& A);

    // Solve given a right-hand side
    virtual void solve(RhsType& x, const RhsType& b) = 0;

    // Storage order
    virtual int order() const = 0;
};

}

}

}

#endif
