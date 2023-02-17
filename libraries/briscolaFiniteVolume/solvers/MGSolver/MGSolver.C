#include "MGSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<>
const char* NamedEnum<MGCycleType,3>::names[] =
{
    "V",
    "W",
    "F"
};

template<class SType, class Type, class MeshType>
void MGSolver<SType,Type,MeshType>::cycle
(
    linearSystem<SType,Type,MeshType>& xEqn,
    meshField<Type,MeshType>& r,
    labelList& visits,
    const label l,
    const label l0,
    const labelList& converged,
    const label nSweepsPre,
    const label nSweepsPost
)
{
    meshField<Type, MeshType>& x = xEqn.x();
    meshField<Type, MeshType>& b = xEqn.b();

    const bool defect = l > l0;

    // Initialize defects to zero

    if (defect)
        x[l] = Zero;

    // Pre-smooth

    this->smooth(xEqn, l, defect, nSweepsPre, converged, omega_);

    // If we are not on the coarsest level, continue to traverse levels.
    // Otherwise, just smooth (could be better to use a direct solver here)

    if (l < x.size()-1)
    {
        for (label rep = 0; rep < levelReps(l,l0,visits); rep++)
        {
            // Don't compute the residual for the finest level during the first
            // repetition

            if (rep > 0 || defect)
                forAll(x[l], d)
                    if (!converged[d])
                        xEqn.residual(r[l][d]);

            // Restrict the current level residual to coarse level

            forAll(x[l], d)
                if (!converged[d])
                    reScheme_->restrict(b[l+1][d], r[l][d], true);

            // Solve the coarse level defect equation

            cycle
            (
                xEqn,
                r,
                visits,
                l+1,
                l0,
                converged,
                nSweepsPre,
                nSweepsPost
            );

            // Prolong the coarse level defect solution and add the correction
            // directly to the current level

            forAll(x[l], d)
                if (!converged[d])
                    proScheme_->prolong(x[l][d], x[l+1][d], plusEqOp<Type>());

            x[l].correctBoundaryConditions(defect);

            // Post-smooth

            this->smooth
            (
                xEqn,
                l,
                defect,
                nSweepsPost,
                converged,
                omega_
            );
        }
    }
    else
    {
        this->smooth
        (
            xEqn,
            l,
            defect,
            Foam::max(nSweepsPost,2),
            converged,
            omega_
        );
    }

    visits[l]++;
}

template<class SType, class Type, class MeshType>
void MGSolver<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const label l,
    const scalar relTol,
    const scalar absTol,
    const label minIter,
    const label maxIter,
    const label nSweepsPre,
    const label nSweepsPost
)
{
    meshField<Type, MeshType>& x = xEqn.x();

    // Correct the system's boundary conditions

    for (label i = l; i < x.size(); i++)
    {
        x[i].correctBoundaryConditions(i > l);
    }

    // Residual field

    meshField<Type, MeshType> r
    (
        IOobject::groupName(x.name(), "residual"),
        x.fvMsh()
    );

    // Residual normalization factors

    const List<Type> normFactors(this->normFactors(xEqn,l));

    // Initial residual

    xEqn.residual(r[l]);

    const List<Type> initialResiduals =
        cmptDivide(gSum(cmptMag(r[l])), normFactors);

    List<Type> currentResiduals(initialResiduals);

    labelList converged =
        this->checkConvergence
        (
            currentResiduals,
            initialResiduals,
            relTol,
            absTol
        );

    label iter = 0;

    labelList visits(x.size());

    while
    (
        (iter < maxIter && Foam::min(converged) == 0)
     || iter < minIter
    )
    {
        // Make sure that we sweep all directions for at least minIter
        // iterations

        if (iter < minIter)
        {
            converged = 0;
        }

        // Cycle

        visits = 0;

        cycle
        (
            xEqn,
            r,
            visits,
            l,
            l,
            converged,
            nSweepsPre,
            nSweepsPost
        );

        // Recompute the residual

        forAll(x[l], d)
            if (!converged[d])
                xEqn.residual(r[l][d]);

        currentResiduals =
            cmptDivide(gSum(cmptMag(r[l])), normFactors);

        converged =
            this->checkConvergence
            (
                currentResiduals,
                initialResiduals,
                relTol,
                absTol
            );

        iter++;
    }

    // Print solver statistics

    this->printSolverStats
    (
        MGSolver<SType,Type,MeshType>::typeName,
        x.name(),
        initialResiduals,
        currentResiduals,
        iter
    );
}

template<class SType, class Type, class MeshType>
MGSolver<SType,Type,MeshType>::MGSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>(dict,fvMsh),
    smoother_
    (
        smootherTypeNames
        [
            dict.lookupOrDefault<word>
            (
                "smoother",
                Pstream::parRun()
              ? smootherTypeNames[LEXGS]
              : smootherTypeNames[RBGS]
            )
        ]
    ),
    omega_
    (
        (smoother_ == RBGS || smoother_ == LEXGS)
      ? dict.lookupOrDefault<scalar>("omega", 1.0)
      : dict.lookupOrDefault<scalar>("omega", 0.8)
    ),
    nSweepsPre_(dict.lookupOrDefault<label>("nSweepsPre", 0)),
    nSweepsPost_(dict.lookupOrDefault<label>("nSweepsPost", 2)),
    cycleType_
    (
        MGCycleTypeNames[dict.lookupOrDefault<word>("cycleType", "F")]
    ),
    proScheme_
    (
        prolongationScheme<Type,MeshType>::New
        (
            dict.lookupOrDefault<word>("prolong", "linear"),
            fvMsh
        )
    ),
    reScheme_
    (
        restrictionScheme<Type,MeshType>::New
        (
            dict.lookupOrDefault<word>("restrict", "linear"),
            fvMsh
        )
    )
{
    // Set smoother function pointer

    if (smoother_ == RBGS)
    {
        smooth = &this->RBGS;
    }
    else if (smoother_ == LEXGS)
    {
        smooth = &this->LEXGS;
    }
    else
    {
        smooth = &this->JAC;
    }
}

template<class SType, class Type, class MeshType>
void MGSolver<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn
)
{
    meshField<Type, MeshType>& x = xEqn.x();

    // Solve starting from level 0

    this->solve(xEqn,0);

    // Restrict solution from fine to coarse

    reScheme_->restrict(x);
}

}

}

}
