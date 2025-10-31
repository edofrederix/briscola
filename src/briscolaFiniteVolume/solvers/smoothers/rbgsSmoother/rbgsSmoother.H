#ifndef rbgsSmoother_H
#define rbgsSmoother_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class rbgsSmoother
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    using solver<SType,Type,MeshType>::smoother::G;
    using solver<SType,Type,MeshType>::smoother::N;

    // Static sweep function on a direction
    static void Sweep
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label d
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
    TypeName("rbgs");

    // Construct from dictionary and mesh
    rbgsSmoother(const dictionary&, const fvMesh&);

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
    #include "rbgsSmoother.C"
#endif

#endif
