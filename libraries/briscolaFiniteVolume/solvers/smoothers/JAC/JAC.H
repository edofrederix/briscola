#ifndef JAC_H
#define JAC_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class JAC
:
    public solver<SType,Type,MeshType>::smoother
{
    // Relaxation factor
    const scalar omega_;

public:

    // Runtime type information
    TypeName("JAC");

    // Construct from dictionary and mesh
    JAC(const dictionary&, const fvMesh&);

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
    #include "JAC.C"
#endif

#endif
