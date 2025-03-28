#ifndef splitPoissonSolver_H
#define splitPoissonSolver_H

#include "PoissonSolver.H"
#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class splitPoissonSolver
:
    public PoissonSolver<SType,Type,MeshType>
{
protected:

    // An optional pre-conditioner to improve the solution extrapolation
    autoPtr<PoissonSolver<SType,Type,MeshType>> preSolverPtr_;

    // The actual Poisson solver
    autoPtr<PoissonSolver<SType,Type,MeshType>> solverPtr_;

    // Cached modified right-hand side
    autoPtr<meshField<Type,MeshType>> bHatPtr_;

    // Number of corrections (defaults to 1)
    label nCorr_;

    // Use a solution extrapolation as initial prediction (defaults to true)
    Switch prediction_;

public:

    // Runtime type information
    TypeName("split");

    // Construct from dictionary and mesh
    splitPoissonSolver
    (
        const dictionary&,
        const fvMesh&
    );

    // Enable the computation of the flux
    virtual void enableFlux()
    {
        solverPtr_->enableFlux();
        PoissonSolver<SType,Type,MeshType>::enableFlux();
    }

    // Disable the computation of the flux
    virtual void disableFlux()
    {
        solverPtr_->disableFlux();
        PoissonSolver<SType,Type,MeshType>::disableFlux();
    }

    // Solve function
    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const meshField<faceScalar,MeshType>* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    );
};

}

}

}

#ifdef NoRepository
    #include "splitPoissonSolver.C"
#endif

#endif
