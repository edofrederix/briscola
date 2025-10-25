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

    using solver<SType,Type,MeshType>::smoother::G;
    using solver<SType,Type,MeshType>::smoother::N;

    // Static sweep function on a direction
    static inline void Sweep
    (
        meshDirection<Type,MeshType>& x,
        const meshDirection<SType,MeshType>& A,
        const meshDirection<Type,MeshType>& b,
        const meshDirection<label,MeshType>& f
    );

    // Static smooth function
    static void Smooth
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label sweeps,
        const labelList& converged
    );

    // Runtime type information
    TypeName("RBGS");

    // Construct from dictionary and mesh
    RBGS(const dictionary&, const fvMesh&);

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
    #include "RBGS.C"
#endif

#endif
