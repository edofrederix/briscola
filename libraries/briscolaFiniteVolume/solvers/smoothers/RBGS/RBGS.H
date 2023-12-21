#ifndef RBGS_H
#define RBGS_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class RBGS
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    // Runtime type information
    TypeName("RBGS");

    // Construct from dictionary and mesh
    RBGS(const dictionary&, const fvMesh&);

    // Smooth at given level
    virtual void smooth
    (
        linearSystem<SType,Type,MeshType>& sys,
        List<Type>& xi,
        const label l,
        const label sweeps,
        const labelList& converged
    );
};

}

}

}

#ifdef NoRepository
    #include "RBGS.C"
#endif

#endif
