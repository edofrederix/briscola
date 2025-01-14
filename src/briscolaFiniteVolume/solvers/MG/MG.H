#ifndef MG_H
#define MG_H

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

// Enumeration of coarse level modes
enum MGCoarseMode
{
    SMOOTH,
    DIRECT
};

// Cycle type names
static const NamedEnum<MGCycleType,3> MGCycleTypeNames;

// Coarse level modes
static const NamedEnum<MGCoarseMode,2> MGCoarseModeNames;

template<class SType, class Type, class MeshType>
class MG
:
    public solver<SType,Type,MeshType>
{
protected:

    // Number of pre-smoother sweeps. Defaults to 0
    label nSweepsPre_;

    // Number of post-smoother sweeps. Defaults to 2
    label nSweepsPost_;

    // Cycle type
    MGCycleType cycleType_;

    // Coarse level mode
    MGCoarseMode coarseMode_;

    // Coarsest level, counted in reversed order with zero being the coarsest
    label coarseLevel_;

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
        linearSystem<SType,Type,MeshType>& sys,
        const scalar relTol,
        const scalar absTol,
        const label minIter,
        const label maxIter,
        const label nSweepsPre,
        const label nSweepsPost,
        const label coarseLevel
    );

    // Multigrid cycle for given mesh level
    virtual void cycle
    (
        linearSystem<SType,Type,MeshType>& sys,
        meshField<Type,MeshType>& r,
        labelList& visits,
        const label l,
        const labelList& converged,
        const label nSweepsPre,
        const label nSweepsPost,
        const label coarseLevel
    );

public:

    // Runtime type information
    TypeName("MG");

    // Construct from dictionary and mesh
    MG(const dictionary&, const fvMesh&);

    // Solve the given system with optional constant matrix indicator
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>& sys,
        const bool constMatrix = false
    );

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
    #include "MG.C"
#endif

#endif
