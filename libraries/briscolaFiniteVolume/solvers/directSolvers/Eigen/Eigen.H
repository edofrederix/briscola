#ifndef Eigen_H
#define Eigen_H

#include "solver.H"
#include "EigenSolver.H"
#include "EigenLinearSystem.H"

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

    // Eigen linear system
    autoPtr<EigenLinearSystem<SType,Type,MeshType>> APtr_;

    // Solvers. By default the SparseLU solver is used.
    PtrList<EigenSolverBase> solverPtrs_;

public:

    // Runtime type information
    TypeName("Eigen");

    // Construct from dictionary and mesh
    Eigen(const dictionary&, const fvMesh&, const label = 0);

    // Prepare
    virtual void prepare(linearSystem<SType,Type,MeshType>&);

    // Solve
    virtual void solve(linearSystem<SType,Type,MeshType>&);
};

}

}

}

#ifdef NoRepository
    #include "Eigen.C"
#endif

#endif
