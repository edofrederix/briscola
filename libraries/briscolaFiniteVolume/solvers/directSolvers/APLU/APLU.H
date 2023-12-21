#ifndef APLU_H
#define APLU_H

#include "solver.H"

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

    // Processor mesh shapes and sizes, per direction. Only fully set on master.
    List<List<labelVector>> shapes_;
    List<List<label>> sizes_;

    // Total data size of communication to collect all cells in all directions
    // into one send/recv
    labelList commSizes_;

    // Stencil indices, for each direction, processor and cell
    List<List<List<labelList>>> indices_;

    // Linear system sizes
    List<label> n_;

public:

    // Runtime type information
    TypeName("APLU");

    // Construct from dictionary and mesh
    APLU(const dictionary&, const fvMesh&, const label = 0);

    // Return the indices
    const List<List<List<labelList>>>& indices() const { return indices_; }

    // Solve the given system, indicating of singular system augmentation should
    // be applied
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
    #include "APLU.C"
#endif

#endif
