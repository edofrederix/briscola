#ifndef iterative_H
#define iterative_H

#include "MG.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

/**

@brief Simple iterative solver using smoothers

Relies on MG to provide the smoother and read parameters. This solver is only
useful for diagonally dominant matrices, e.g., for the Navier-Stokes equations
at large Re.

*/

template<class SType, class Type, class MeshType>
class iterative
:
    public MG<SType,Type,MeshType>
{
protected:

    // Number of smoother sweeps in each smoother call before checking the
    // residual (defaults to 1)
    label nSweeps_;

    // Solve system and provide parameters
    void solve
    (
        linearSystem<SType,Type,MeshType>& sys,
        const scalar relTol,
        const scalar absTol,
        const label minIter,
        const label maxIter,
        const label nSweeps
    );

public:

    // Runtime type information
    TypeName("iterative");

    // Construct from dictionary and mesh
    iterative(const dictionary&, const fvMesh&);

    // Solve the given system with optional constant matrix indicator
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>& sys,
        const bool constMatrix = false
    );
};

}

}

}

#ifdef NoRepository
    #include "iterative.C"
#endif

#endif
