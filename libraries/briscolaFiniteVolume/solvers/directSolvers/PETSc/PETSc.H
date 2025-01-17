#ifndef PETSc_H
#define PETSc_H

#include "solver.H"
#include "PETScSolverBase.H"
#include "PETScLinearSystem.H"
#include <petscmat.h>

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class PETSc
:
    public solver<SType,Type,MeshType>::directSolver
{
protected:

    // PETSc linear system
    autoPtr<PETScLinearSystem<SType,Type,MeshType>> APtr_;

    // Solvers
    PtrList<PETScSolverBase> solverPtrs_;

    // Number of aggregation parts (defaults to 1 and forced to 1 if the
    // selected PETSc solver is not parallel)
    label nAggregationParts_;

    // Print solver statistics (defaults to false)
    Switch printStats_;

public:

    // Runtime type information
    TypeName("PETSc");

    // Construct from dictionary and mesh
    PETSc(const dictionary&, const fvMesh&, const label = 0);

    // Destructor
    ~PETSc();

    // Prepare
    virtual void prepare(linearSystem<SType,Type,MeshType>&);

    // Solve
    virtual void solve(linearSystem<SType,Type,MeshType>&);
};

}

}

}

#ifdef NoRepository
    #include "PETSc.C"
#endif

#endif
