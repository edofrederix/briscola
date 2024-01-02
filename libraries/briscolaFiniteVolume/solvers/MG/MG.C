#include "MG.H"

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

template<>
const char* NamedEnum<MGCoarseMode,2>::names[] =
{
    "smooth",
    "direct"
};

template<class SType, class Type, class MeshType>
void MG<SType,Type,MeshType>::cycle
(
    linearSystem<SType,Type,MeshType>& sys,
    meshField<Type,MeshType>& r,
    labelList& visits,
    const label l,
    const List<bool>& singular,
    const labelList& converged,
    const label nSweepsPre,
    const label nSweepsPost,
    const label coarseLevel
)
{
    meshField<Type, MeshType>& x = sys.x();
    meshField<Type, MeshType>& b = sys.b();

    // Initialize defects to zero

    if (l > 0)
        x[l] = Zero;

    // Only apply singular system augmentation at the coarsest level and if any
    // of the directions is actually singular. Initialize the auxiliary unknown
    // to zero.

    List<Type> xi(0);

    if (l == x.size()-1 && sum(singular))
    {
        xi.setSize(x[l].size(), Zero);
    }

    // Pre-smooth

    this->smoothPtr_->smooth
    (
        sys,
        xi,
        l,
        nSweepsPre,
        converged
    );

    // If we are not on the coarsest level, continue to traverse levels.
    // Otherwise, just smooth (could be better to use a direct solver here)

    if (l < x.size()-1-coarseLevel)
    {
        for (label rep = 0; rep < levelReps(l,visits); rep++)
        {
            // Don't compute the residual for the finest level during the first
            // repetition

            if (rep > 0 || l > 0)
                forAll(x[l], d)
                    if (!converged[d])
                        sys.residual(r[l][d]);

            // Restrict the current level residual to coarse level

            forAll(x[l], d)
                if (!converged[d])
                    reScheme_->restrict(b[l+1][d], r[l][d], true);

            // Solve the coarse level defect equation

            cycle
            (
                sys,
                r,
                visits,
                l+1,
                singular,
                converged,
                nSweepsPre,
                nSweepsPost,
                coarseLevel
            );

            // Prolong the coarse level defect solution and add the correction
            // directly to the current level

            forAll(x[l], d)
                if (!converged[d])
                    proScheme_->prolong(x[l][d], x[l+1][d], plusEqOp<Type>());

            x[l].correctBoundaryConditions();

            // Post-smooth

            this->smoothPtr_->smooth
            (
                sys,
                xi,
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
            this->directSolvePtr_->solve(sys,singular);
        }
        else
        {
            this->smoothPtr_->smooth
            (
                sys,
                xi,
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
    const List<bool>& singular,
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

    // Correct the boundary conditions

    x[0].correctBoundaryConditions();

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

    const List<Type> normFactors(this->normFactors(sys,0));

    // Initial residual

    sys.residual(r[0]);

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
            sys,
            r,
            visits,
            0,
            singular,
            converged,
            nSweepsPre,
            nSweepsPost,
            coarseLevel
        );

        // Recompute the residual

        forAll(x[0], d)
            if (!converged[d])
                sys.residual(r[0][d]);

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
        MGCoarseModeNames[dict.lookupOrDefault<word>("coarseMode", "direct")]
    ),
    coarseLevel_(dict.lookupOrDefault<label>("coarseLevel", 1)),
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
    if (coarseLevel_ > fvMsh.size()-2)
        FatalErrorInFunction
            << "Coarse level is set to high. Should be at most "
            << fvMsh.size()-2 << endl << abort(FatalError);

    this->smoothPtr_.reset
    (
        solver<SType,Type,MeshType>::smoother::New
        (
            this->dict_.template lookupOrDefault<word>
            (
                "smoother",
                Pstream::parRun()
              ? "symmLEXGS"
              : "symmRBGS"
            ),
            this->dict_,
            fvMsh
        ).ptr()
    );

    // Set the coarse level solver

    if (coarseMode_ == DIRECT)
        this->directSolvePtr_.reset
        (
            solver<SType,Type,MeshType>::directSolver::New
            (
                this->dict_.template lookupOrDefault<word>
                (
                    "directSolverType",
                    "Eigen"
                ),
                this->dict_,
                fvMsh,
                fvMsh.msh().size()-1-coarseLevel_
            ).ptr()
        );
}

template<class SType, class Type, class MeshType>
void MG<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
    sys.eliminateGhosts();

    sys.x().makeDeep();
    sys.b().makeDeep();

    List<bool> singular(sys.singular());

    if
    (
        coarseMode_ == DIRECT
     && (
            !constMatrix
         || !this->directSolvePtr_->prepared()
        )
    )
    {
        this->directSolvePtr_->prepare(sys,singular);
    }

    this->solve
    (
        sys,
        singular,
        this->relTol_,
        this->tolerance_,
        this->minIter_,
        this->maxIter_,
        this->nSweepsPre_,
        this->nSweepsPost_,
        this->coarseLevel_
    );

    sys.x().makeShallow();
    sys.b().makeShallow();
}

}

}

}
