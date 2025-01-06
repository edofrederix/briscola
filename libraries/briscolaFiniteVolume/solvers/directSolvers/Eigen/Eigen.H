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

    // Solvers. By default the iterative BiCGSTAB solver is used.
    PtrList<EigenSolverBase> solverPtrs_;

    // Number of aggregation parts (defaults to the number of processors)
    label nAggregationParts_;

    // Minimum number of iterations (defaults to 1)
    label minIter_;

    // Maximum number of iterations (forced to 1 when there's only one
    // aggregation part and defaulting to nAggregationParts otherwise)
    label maxIter_;

    // Relative tolerance (defaults to 1e-3)
    scalar relTol_;

    // Absolute tolerance (defaults to 0)
    scalar tolerance_;

    // Print solver statistics (defaults to false)
    Switch printStats_;

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
