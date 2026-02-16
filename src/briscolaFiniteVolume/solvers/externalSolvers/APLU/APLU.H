#ifndef APLU_H
#define APLU_H

#include "solver.H"
#include "scalarMatrices.H"
#include "linearSystemAggregation.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class APLU
:
    public solver<SType,Type,MeshType>::externalSolver
{
protected:

    // Aggregated linear system
    autoPtr<linearSystemAggregation<SType,Type,MeshType>> APtr_;

    // Dense matrices
    PtrList<scalarSquareMatrix> DPtrs_;

    // Pivot indices
    PtrList<labelList> pivotPtrs_;

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
