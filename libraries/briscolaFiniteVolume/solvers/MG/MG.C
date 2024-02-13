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
    const bool constMatrix,
    const labelList& converged,
    const label nSweepsPre,
    const label nSweepsPost,
    const label coarseLevel
)
{
    meshLevel<Type, MeshType>& xl = sys.x()[l];
    meshLevel<Type, MeshType>& rl = r[l];

    const fvMesh& fvMsh = this->fvMsh_;

    const label nLevels = this->fvMsh_.size();
    const label nDirs = xl.size();

    // Initialize defects to zero

    if (l > 0)
        xl = Zero;

    // Only apply singular system augmentation at the coarsest level and if any
    // of the directions is actually singular. Initialize the auxiliary unknown
    // to zero.

    List<Type> xi(0);

    if (l == nLevels-1 && sum(singular))
        xi.setSize(nDirs, Zero);

    // Pre-smooth

    if (nSweepsPre > 0)
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

    if (l < nLevels-1-coarseLevel)
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

                        if
                        (
                               (fvMsh.immersedBoundaryPresent())
                            && (sys.x().IBC().Jac())
                        )
                        {
                            forAllCells(r[l][d],i,j,k)
                            {
                                if
                                (
                                    fvMsh.IB<MeshType>()
                                        .ghostMask()[l][d](i,j,k)
                                )
                                {
                                    r[l][d](i,j,k) = Zero;
                                }
                            }
                        }
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
                singular,
                constMatrix,
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
            xl.correctImmersedBoundaryConditions();

            // Post-smooth

            if (nSweepsPost > 0)
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
    const bool constMatrix,
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

    const fvMesh& fvMsh = this->fvMsh_;

    // Correct the boundary conditions

    x[0].correctBoundaryConditions();
    x[0].correctImmersedBoundaryConditions();

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

    // Initial residual

    sys.residual(r[0]);

    if
    (
           (fvMsh.immersedBoundaryPresent())
        && (sys.x().IBC().Jac())
    )
    {
        forAllCells(r[0],d,i,j,k)
        {
            if (fvMsh.IB<MeshType>().ghostMask()[0][d](i,j,k))
            {
                r[0][d](i,j,k) = Zero;
            }
        }
    }


    // Residual normalization factors

    const List<Type> normFactors(this->normFactors(sys,r,0));

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
            constMatrix,
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

                if
                (
                       (fvMsh.immersedBoundaryPresent())
                    && (sys.x().IBC().Jac())
                )
                {
                    forAllCells(r[0][d],i,j,k)
                    {
                        if (fvMsh.IB<MeshType>().ghostMask()[0][d](i,j,k))
                        {
                            r[0][d](i,j,k) = Zero;
                        }
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
    coarseLevel_(dict.lookupOrDefault<label>("coarseLevel", 0)),
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
            << "Coarse level is set too high. Should be at most "
            << fvMsh.size()-2 << endl << abort(FatalError);

    this->smoothPtr_.reset
    (
        solver<SType,Type,MeshType>::smoother::New
        (
            this->dict_.template lookupOrDefault<word>
            (
                "smoother",
                "RBGS"
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
                    "directSolver",
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
        constMatrix,
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
