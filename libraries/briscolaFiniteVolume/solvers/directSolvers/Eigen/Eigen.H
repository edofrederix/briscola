#ifndef Eigen_H
#define Eigen_H

#include "solver.H"
#include "globalCellNumbers.H"
#include "globalRowData.H"

#include <eigen3/Eigen/Sparse>

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class Eigen
:
    public solver<SType,Type,MeshType>::directSolver
{
protected:

    // Eigen types
    typedef ::Eigen::SparseMatrix<double> matrixType;
    typedef ::Eigen::MatrixXd rhsType;
    typedef ::Eigen::Triplet<double> tType;
    typedef ::Eigen::SparseLU<matrixType> gSolverType;
    typedef ::Eigen::SimplicialLDLT<matrixType> sSolverType;

    // Sparse matrix
    PtrList<matrixType> APtrs_;

    // Solver
    PtrList<gSolverType> gSolverPtrs_;
    PtrList<sSolverType> sSolverPtrs_;

    // Symmetric system
    List<bool> symmetric_;

    // Global stencil cell numbers
    globalCellNumbers<SType,MeshType> gCellNumbers_;

    // Check if the given sparse matrix is symmetrc
    bool isSymmetric(const matrixType&) const;

public:

    // Runtime type information
    TypeName("Eigen");

    // Construct from dictionary and mesh
    Eigen(const dictionary&, const fvMesh&, const label = 0);

    // Prepare
    virtual void prepare
    (
        const linearSystem<SType,Type,MeshType>&,
        const List<bool>&
    );

    // Solve
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>&,
        const List<bool>&
    );
};

}

}

}

#ifdef NoRepository
    #include "Eigen.C"
#endif

#endif
