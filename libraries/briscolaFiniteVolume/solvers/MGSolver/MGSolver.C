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
    const labelList& converged,
    const label nSweepsPre,
    const label nSweepsPost
)
{
    meshField<Type, MeshType>& x = xEqn.x();
    meshField<Type, MeshType>& b = xEqn.b();

    // Initialize defects to zero

    if (l > 0)
        x[l] = Zero;

    // Pre-smooth

    this->smooth(xEqn, l, nSweepsPre, converged, omega_, this->IB_);

    // If we are not on the coarsest level, continue to traverse levels.
    // Otherwise, just smooth (could be better to use a direct solver here)

    if (l < x.size()-1)
    {
        for (label rep = 0; rep < levelReps(l,visits); rep++)
        {
            // Don't compute the residual for the finest level during the first
            // repetition

            if (rep > 0 || l > 0)
                forAll(x[l], d)
                    if (!converged[d])
                    {
                        xEqn.residual(r[l][d]);
                        if (this->IB_ != nullptr)
                        {
                            forAllCells(r[l][d],i,j,k)
                            {
                                if (this->IB_->ghostMask()[l][d](i,j,k) == 1)
                                {
                                    r[l][d](i,j,k) = Zero;
                                }
                            }
                        }
                    }

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
                converged,
                nSweepsPre,
                nSweepsPost
            );

            // Prolong the coarse level defect solution and add the correction
            // directly to the current level

            forAll(x[l], d)
                if (!converged[d])
                    proScheme_->prolong(x[l][d], x[l+1][d], plusEqOp<Type>());

            x[l].correctCommBoundaryConditions();
            this->correctImmersedBoundaryConditions(this->IB_,xEqn,l);

            // Post-smooth

            this->smooth
            (
                xEqn,
                l,
                nSweepsPost,
                converged,
                omega_,
                this->IB_
            );
        }
    }
    else
    {
        this->smooth
        (
            xEqn,
            l,
            Foam::max(nSweepsPost,2),
            converged,
            omega_,
            this->IB_
        );
    }

    visits[l]++;
}

template<class SType, class Type, class MeshType>
void MGSolver<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const scalar relTol,
    const scalar absTol,
    const label minIter,
    const label maxIter,
    const label nSweepsPre,
    const label nSweepsPost
)
{
    meshField<Type, MeshType>& x = xEqn.x();

    // Correct the boundary conditions

    x[0].correctCommBoundaryConditions();
    this->correctImmersedBoundaryConditions(this->IB_,xEqn,0);

    // Residual field

    meshField<Type, MeshType> r
    (
        IOobject::groupName(x.name(), "residual"),
        x.fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    );

    // Residual normalization factors

    const List<Type> normFactors(this->normFactors(xEqn,0));

    // Initial residual

    xEqn.residual(r[0]);

    if (this->IB_ != nullptr)
    {
        forAllDirections(r[0],d,i,j,k)
        {
            if (this->IB_->ghostMask()[0][d](i,j,k) == 1)
            {
                r[0][d](i,j,k) = Zero;
            }
        }
    }

    const List<Type> initialResiduals =
        cmptDivide(gSum(cmptMag(r[0])), normFactors);

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
            0,
            converged,
            nSweepsPre,
            nSweepsPost
        );

        // Recompute the residual

        forAll(x[0], d)
            if (!converged[d])
            {
                xEqn.residual(r[0][d]);
                if (this->IB_ != nullptr)
                {
                    forAllCells(r[0][d],i,j,k)
                    {
                        if (this->IB_->ghostMask()[0][d](i,j,k) == 1)
                        {
                            r[0][d](i,j,k) = Zero;
                        }
                    }
                }
            }

        currentResiduals =
            cmptDivide(gSum(cmptMag(r[0])), normFactors);

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
        this->typeName,
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
    const fvMesh& fvMsh,
    immersedBoundary<Type,MeshType>* IB
)
:
    solver<SType,Type,MeshType>(dict,fvMsh,IB),
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
    xEqn.eliminateGhosts();

    xEqn.x().makeDeep();
    xEqn.b().makeDeep();

    this->solve
    (
        xEqn,
        this->relTol_,
        this->tolerance_,
        this->minIter_,
        this->maxIter_,
        this->nSweepsPre_,
        this->nSweepsPost_
    );

    xEqn.x().makeShallow();
    xEqn.b().makeShallow();
}

}

}

}
