#ifndef symmRBGS_H
#define symmRBGS_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class symmRBGS
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    // Runtime type information
    TypeName("symmRBGS");

    // Construct from dictionary and mesh
    symmRBGS(const dictionary&, const fvMesh&);

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
    #include "symmRBGS.C"
#endif

#endif
