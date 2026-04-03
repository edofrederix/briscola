#include "iterative.H"
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
void iterative<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const scalar relTol,
    const scalar absTol,
    const label minIter,
    const label maxIter,
    const label nSweeps
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

    if (sigFpeEnabled())
        r = Zero;

    // Initial residual without boundary correction

    sys.template residual<false>(r[0]);

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

        // Smooth

        this->Smooth(sys, 0, nSweeps, converged);

        // Recompute the residual

        sys.residual(r[0]);

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
iterative<SType,Type,MeshType>::iterative
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    MG<SType,Type,MeshType>(dict,fvMsh),
    nSweeps_(dict.lookupOrDefault<label>("nSweeps", 1))
{}

template<class SType, class Type, class MeshType>
void iterative<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
    if (SType::nCsComponents > 1)
        sys.eliminateGhosts();

    sys.setForcingMask();

    if (sys.diagonal())
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            diagonalSmoother<SType,Type,MeshType>::Sweep(sys, 0, d);

        sys.x()[0].correctBoundaryConditions();

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
        this->setSingularityConstraint(sys, 0);

        this->solve
        (
            sys,
            this->relTol_,
            this->tolerance_,
            this->minIter_,
            this->maxIter_,
            this->nSweeps_
        );

        // Correct eliminated ghosts

        sys.x()[0].correctEliminatedBoundaryConditions();
    }
}

}

}

}
