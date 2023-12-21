#ifndef symmLEXGS_H
#define symmLEXGS_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class symmLEXGS
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    // Runtime type information
    TypeName("symmLEXGS");

    // Construct from dictionary and mesh
    symmLEXGS(const dictionary&, const fvMesh&);

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
    #include "symmLEXGS.C"
#endif

#endif
