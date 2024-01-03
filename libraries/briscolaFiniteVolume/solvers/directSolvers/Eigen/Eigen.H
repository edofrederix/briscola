#ifndef Eigen_H
#define Eigen_H

#include "solver.H"
#include "globalCellNumbers.H"
#include "globalRowData.H"
#include "EigenSolver.H"

#ifdef MKL
#define EIGEN_USE_MKL_ALL
#endif

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

    // Global stencil cell numbers
    globalCellNumbers<SType,MeshType> gCellNumbers_;

    // Sparse matrices
    PtrList<EigenSolver::matrixType> APtrs_;

    // Solvers
    PtrList<EigenSolver> solverPtrs_;

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
