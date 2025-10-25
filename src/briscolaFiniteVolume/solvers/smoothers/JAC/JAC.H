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

    // Relaxation factor
    static const scalar omega_;

    // Static smooth function
    static void Smooth
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label sweeps,
        const labelList& converged
    );

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
