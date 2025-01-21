#ifndef LEXGS_H
#define LEXGS_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class LEXGS
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    // Static smooth function
    static void Smooth
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label sweeps,
        const labelList& converged
    );

    // Runtime type information
    TypeName("LEXGS");

    // Construct from dictionary and mesh
    LEXGS(const dictionary&, const fvMesh&);

    // Smooth at given level
    virtual void smooth
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label sweeps,
        const labelList& converged
    );
};

}

}

}

#ifdef NoRepository
    #include "LEXGS.C"
#endif

#endif
