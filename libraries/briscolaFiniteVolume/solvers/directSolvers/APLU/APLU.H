#ifndef APLU_H
#define APLU_H

#include "solver.H"
#include "globalCellNumbers.H"
#include "globalRowData.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Aggregated Pivoted LU solver

template<class SType, class Type, class MeshType>
class APLU
:
    public solver<SType,Type,MeshType>::directSolver
{
protected:

    // Dense matrix
    PtrList<scalarSquareMatrix> APtrs_;

    // Pivot indices
    PtrList<labelList> pivotPtrs_;

    // Global stencil cell numbers
    globalCellNumbers<SType,MeshType> gCellNumbers_;

public:

    // Runtime type information
    TypeName("APLU");

    // Construct from dictionary and mesh
    APLU(const dictionary&, const fvMesh&, const label = 0);

    // Prepare
    virtual void prepare(linearSystem<SType,Type,MeshType>&);

    // Solve
    virtual void solve(linearSystem<SType,Type,MeshType>&);
};

}

}

}

#ifdef NoRepository
    #include "APLU.C"
#endif

#endif
