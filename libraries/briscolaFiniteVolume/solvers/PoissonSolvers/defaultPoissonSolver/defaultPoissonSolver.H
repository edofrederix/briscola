#ifndef defaultPoissonSolver_H
#define defaultPoissonSolver_H

#include "PoissonSolver.H"
#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class defaultPoissonSolver
:
    public PoissonSolver<SType,Type,MeshType>
{
protected:

    // The actual Poisson solver
    autoPtr<solver<SType,Type,MeshType>> solverPtr_;

    // Cached linear system pointer
    autoPtr<linearSystem<SType,Type,MeshType>> sysPtr_;

public:

    // Runtime type information
    TypeName("default");

    // Construct from dictionary and mesh
    defaultPoissonSolver
    (
        const word,
        const dictionary&,
        const fvMesh&
    );

    // Solve function
    virtual void solve
    (
        meshField<Type,MeshType>&,
        const meshField<Type,MeshType>*,
        const meshField<lowerFaceScalar,MeshType>*,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "defaultPoissonSolver.C"
#endif

#endif
