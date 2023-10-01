#ifndef MGSolver_H
#define MGSolver_H

#include "solver.H"
#include "prolongationScheme.H"
#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Enumeration of cycle types
enum MGCycleType
{
    VCYCLE,
    WCYCLE,
    FCYCLE
};

// Cycle type names
static const NamedEnum<MGCycleType,3> MGCycleTypeNames;

template<class SType, class Type, class MeshType>
class MGSolver
:
    public solver<SType,Type,MeshType>
{
protected:

    // Smoother type
    smootherType smoother_;

    // Smoothing relaxation factor. When not set, it defaults to 1 for GS and
    // 0.8 for JAC smoothers
    scalar omega_;

    // Number of pre-smoother sweeps. Defaults to 0
    label nSweepsPre_;

    // Number of post-smoother sweeps. Defaults to 2
    label nSweepsPost_;

    // Cycle type
    MGCycleType cycleType_;

    // Number of level repetitions. Use 1 for a V-cycle, use 2 for a W-cycle.
    // Defaults to 2
    label levelReps_;

    // Prolongation operator. By default the midPoint scheme is used
    autoPtr<prolongationScheme<Type,MeshType>> proScheme_;

    // Restriction operator. By default the average scheme is used
    autoPtr<restrictionScheme<Type,MeshType>> reScheme_;

    // Solve system and provide parameters
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>& xEqn,
        const scalar relTol,
        const scalar absTol,
        const label minIter,
        const label maxIter,
        const label nSweepsPre,
        const label nSweepsPost
    );

    // Multigrid cycle for given mesh level
    virtual void cycle
    (
        linearSystem<SType,Type,MeshType>& xEqn,
        meshField<Type,MeshType>& r,
        labelList& visits,
        const label l,
        const labelList& converged,
        const label nSweepsPre,
        const label nSweepsPost
    );

    // Smoother function pointer
    void (*smooth)
    (
        linearSystem<SType,Type,MeshType>&,
        const label,
        const label,
        const labelList&,
        const scalar
    );

public:

    // Runtime type information
    TypeName("MG");

    // Construct from dictionary and mesh
    MGSolver(const dictionary&, const fvMesh&);

    // Solve the given system
    virtual void solve(linearSystem<SType,Type,MeshType>&);

    // Number of level repetitions
    inline label levelReps
    (
        const label l,
        const labelList& visits
    )
    {
        if (cycleType_ == VCYCLE)
        {
            return 1;
        }
        else if (cycleType_ == WCYCLE)
        {
            return 1 + (l > 0);
        }
        else
        {
            return 1 + (l > 0 && visits[l] == 0);
        }
    }
};

}

}

}

#ifdef NoRepository
    #include "MGSolver.C"
#endif

#endif
