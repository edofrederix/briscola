#ifndef Krylov_H
#define Krylov_H

#include "solver.H"
#include <petscksp.h>
#include "FastPtrList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class Krylov
:
    public solver<SType,Type,MeshType>
{
protected:

    // PETSc Krylov solver type, see
    // https://petsc.org/main/manualpages/KSP/KSPType
    word solverType_;

    // PETSc matrices
    FastPtrList<Mat> matrices_;

    // PETSc solvers
    FastPtrList<KSP> solvers_;

    // Preparation bool
    bool prepared_;

    // Prepare matrix and solver
    virtual void prepare
    (
        linearSystem<SType,Type,MeshType>& sys
    );

    // Solve system and provide parameters
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>& sys,
        const scalar relTol,
        const scalar absTol,
        const label maxIter
    );

public:

    // Runtime type information
    TypeName("Krylov");

    // Construct from dictionary and mesh
    Krylov(const dictionary&, const fvMesh&);

    // Destructor
    ~Krylov();

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
    #include "Krylov.C"
#endif

#endif
