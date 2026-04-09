#ifndef jacSmoother_H
#define jacSmoother_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class jacSmoother
:
    public solver<SType,Type,MeshType>::smoother
{
public:

    using solver<SType,Type,MeshType>::smoother::G;
    using solver<SType,Type,MeshType>::smoother::N;

    // Relaxation factor
    static const scalar omega_;

    // Static sweep function on a direction with optional padding
    static void Sweep
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label d,
        const faceLabel& paddingMask,
        const label padding = 0
    );

    // Static boundary sweep function on a direction given an active boundary
    // mask and optional padding (defaults to one)
    static void SweepBoundary
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label d,
        const faceLabel& paddingMask,
        const label padding = 1
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
    TypeName("jac");

    // Construct from dictionary and mesh
    jacSmoother(const dictionary&, const fvMesh&);

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
    #include "jacSmoother.C"
#endif

#endif
