#include "MG.H"
#include "diagonal.H"

#include "diagonalSmoother.H"
#include "rbgsSmoother.H"
#include "lexgsSmoother.H"
#include "jacSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void MG<SType,Type,MeshType>::cycle
(
    linearSystem<SType,Type,MeshType>& sys,
    meshField<Type,MeshType>& r,
    labelList& visits,
    const label l,
    const labelList& converged,
    const label nSweepsPre,
    const label nSweepsPost,
    const label coarseLevel
)
{
    meshLevel<Type, MeshType>& xl = sys.x()[l];
    meshLevel<Type, MeshType>& rl = r[l];

    // On the coarsest level force the first global cell to zero for a
    // singular matrix

    List<bool> singular(sys.singular());

    forAll(singular, d)
        if (l == coarseLevel && singular[d] && Pstream::master())
            sys.b()[l][d](0,0,0) = Zero;

    // Initialize defects to zero

    if (l > 0)
        xl = Zero;

    // Pre-smooth

    if (nSweepsPre > 0)
        this->Smooth
        (
            sys,
            l,
            nSweepsPre,
            converged
        );

    // If we are not on the coarsest level, continue to traverse levels.
    // Otherwise solve the coarsest level.

    if (l < coarseLevel)
    {
        meshLevel<Type, MeshType>& xlCoarse = sys.x()[l+1];
        meshLevel<Type, MeshType>& blCoarse = sys.b()[l+1];

        for (label rep = 0; rep < levelReps(l,visits); rep++)
        {
            // Don't compute the residual for the finest level during the first
            // repetition

            if (rep > 0 || l > 0)
            {
                forAll(xl, d)
                {
                    if (!converged[d])
                    {
                        sys.residual(rl[d]);
                    }
                }
            }

            // Restrict the current level residual to coarse level

            forAll(xl, d)
                if (!converged[d])
                    reScheme_->restrict(blCoarse[d], rl[d], true);

            // Solve the coarse level defect equation

            cycle
            (
                sys,
                r,
                visits,
                l+1,
                converged,
                nSweepsPre,
                nSweepsPost,
                coarseLevel
            );

            // Prolong the coarse level defect solution and add the correction
            // directly to the current level

            forAll(xl, d)
                if (!converged[d])
                    proScheme_->prolong(xl[d], xlCoarse[d], plusEqOp<Type>());

            xl.correctBoundaryConditions();

            // Post-smooth

            if (nSweepsPost > 0)
                this->Smooth
                (
                    sys,
                    l,
                    nSweepsPost,
                    converged
                );
        }
    }
    else
    {
        if (coarseMode_ == DIRECT)
        {
            this->coarseSolverPtr_->solve(sys);
        }
        else
        {
            this->Smooth
            (
                sys,
                l,
                Foam::max(nSweepsPost,2),
                converged
            );
        }
    }

    visits[l]++;
}

template<class SType, class Type, class MeshType>
void MG<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const scalar relTol,
    const scalar absTol,
    const label minIter,
    const label maxIter,
    const label nSweepsPre,
    const label nSweepsPost,
    const label coarseLevel
)
{
    meshField<Type, MeshType>& x = sys.x();

    // Residual field

    meshField<Type, MeshType> r
    (
        IOobject::groupName(x.name(), "residual"),
        x.fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    );

    // Initial residual

    sys.residual(r[0]);

    // Residual normalization factors

    const List<Type> normFactors(this->normFactors(sys,r[0]));

    const List<Type> initialResiduals =
        cmptDivide
        (
            Foam::cmptSqrt(gSum(cmptSqr(r[0]))),
            normFactors
        );

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
            sys,
            r,
            visits,
            0,
            converged,
            nSweepsPre,
            nSweepsPost,
            coarseLevel
        );

        // Recompute the residual

        forAll(x[0], d)
        {
            if (!converged[d])
            {
                sys.residual(r[0][d]);
            }
        }

        currentResiduals =
            cmptDivide
            (
                Foam::cmptSqrt(gSum(cmptSqr(r[0]))),
                normFactors
            );

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
MG<SType,Type,MeshType>::MG
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>(dict,fvMsh),
    nSweepsPre_(dict.lookupOrDefault<label>("nSweepsPre", 0)),
    nSweepsPost_(dict.lookupOrDefault<label>("nSweepsPost", 2)),
    cycleType_
    (
        MGCycleTypeNames[dict.lookupOrDefault<word>("cycleType", "F")]
    ),
    coarseMode_
    (
        MGCoarseModeNames[dict.lookupOrDefault<word>("coarseMode", "smooth")]
    ),
    coarseLevel_(dict.lookupOrDefault<label>("coarseLevel", 0)),
    proScheme_
    (
        prolongationScheme<Type,MeshType>::New
        (
            fvMsh,
            dict.lookupOrDefault<word>("prolong", "linear")
        )
    ),
    reScheme_
    (
        restrictionScheme<Type,MeshType>::New
        (
            fvMsh,
            dict.lookupOrDefault<word>("restrict", "linear")
        )
    )
{
    if (coarseLevel_ > fvMsh.size()-2)
        FatalErrorInFunction
            << "Coarse level is set too high. Should be at most "
            << fvMsh.size()-2 << endl << abort(FatalError);

    // Set the smoother of choice only if the stencil is not diagonal

    if (SType::nCsComponents > 1)
        this->smoothPtr_.reset
        (
            solver<SType,Type,MeshType>::smoother::New
            (
                this->dict_.template lookupOrDefault<word>
                (
                    "smoother",
                    "rbgs"
                ),
                this->dict_,
                fvMsh
            ).ptr()
        );

    // Set smoother pointer to avoid virtual function call

    if (this->smoothPtr_.valid())
    {
        if
        (
            this->smoothPtr_->type()
         == rbgsSmoother<SType,Type,MeshType>::typeName
        )
        {
            Smooth =
                &dynamic_cast<rbgsSmoother<SType,Type,MeshType>*>
                (
                    &this->smoothPtr_()
                )->Smooth;
        }
        else if
        (
            this->smoothPtr_->type()
         == lexgsSmoother<SType,Type,MeshType>::typeName
        )
        {
            Smooth =
                &dynamic_cast<lexgsSmoother<SType,Type,MeshType>*>
                (
                    &this->smoothPtr_()
                )->Smooth;
        }
        else if
        (
            this->smoothPtr_->type()
         == jacSmoother<SType,Type,MeshType>::typeName
        )
        {
            Smooth =
                &dynamic_cast<jacSmoother<SType,Type,MeshType>*>
                (
                    &this->smoothPtr_()
                )->Smooth;
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported smoother type" << endl << abort(FatalError);
        }
    }

    // Set the coarse level solver if requested to do so, and only if the
    // stencil is not diagonal

    if (coarseMode_ == DIRECT && SType::nCsComponents > 1)
    {
        if (!this->dict_.found("coarseSolver"))
            this->dict_.add("coarseSolver", dictionary());

        dictionary& subDict = this->dict_.subDict("coarseSolver");

        subDict.lookupOrAddDefault("type", word("PETSc"));
        subDict.lookupOrAddDefault("relTol", scalar(1e-3));

        const label coarseLevel = fvMsh.msh().size() - coarseLevel_ - 1;

        this->coarseSolverPtr_.reset
        (
            solver<SType,Type,MeshType>::externalSolver::New
            (
                this->dict_.subDict("coarseSolver"),
                fvMsh,
                coarseLevel
            ).ptr()
        );
    }
}

template<class SType, class Type, class MeshType>
void MG<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
    if (SType::nCsComponents > 1)
        sys.eliminateGhosts();

    sys.setForcingMask();

    if
    (
        SType::nCsComponents == 1
     || sum(labelList(sys.diagonal())) == MeshType::numberOfDirections
    )
    {
        sys.x().makeShallow();
        sys.b().makeShallow();

        diagonalSmoother<SType,Type,MeshType>::Smooth
        (
            sys,
            0,
            1,
            labelList(MeshType::numberOfDirections, 0)
        );

        this->printSolverStats
        (
            diagonal<SType,Type,MeshType>::typeName,
            sys.x().name(),
            List<Type>(MeshType::numberOfDirections, pTraits<Type>::one),
            List<Type>(MeshType::numberOfDirections, Zero),
            0
        );
    }
    else
    {
        const label coarseLevel = sys.A().size() - this->coarseLevel_ - 1;

        this->setSingularityConstraint(sys, coarseLevel);

        sys.x().makeDeep();
        sys.b().makeDeep();

        if
        (
            coarseMode_ == DIRECT
        && (
                !constMatrix
             || !this->coarseSolverPtr_->prepared()
            )
        )
        {
            this->coarseSolverPtr_->prepare(sys);
        }

        this->solve
        (
            sys,
            this->relTol_,
            this->tolerance_,
            this->minIter_,
            this->maxIter_,
            this->nSweepsPre_,
            this->nSweepsPost_,
            coarseLevel
        );

        sys.x().makeShallow();
        sys.b().makeShallow();
    }
}

}

}

}
