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

    // Correct the boundary conditions

tic(4)
    x[0].correctBoundaryConditions();
toc(4)
    // Residual field
tic(5)

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
toc(5)
tic(6)

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
toc(6)

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

tic(7)
        this->Smooth(sys, 0, nSweeps, converged);
toc(7)

        // Recompute the residual

tic(8)
        forAll(x[0], d)
            if (!converged[d])
                sys.residual(r[0][d]);
toc(8)

tic(9)
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
toc(9)

        iter++;
    }

    // Print solver statistics
tic(10)
    this->printSolverStats
    (
        this->typeName,
        x.name(),
        initialResiduals,
        currentResiduals,
        iter
    );
toc(10)
}

template<class SType, class Type, class MeshType>
iterative<SType,Type,MeshType>::iterative
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    MG<SType,Type,MeshType>(dict,fvMsh),
    initTicTocConstructor(12),
    nSweeps_(dict.lookupOrDefault<label>("nSweeps", 1))
{}

template<class SType, class Type, class MeshType>
void iterative<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
tic(0)
    if (SType::nCsComponents > 1)
        sys.eliminateGhosts();
toc(0)
tic(1)
    sys.x().makeShallow();
    sys.b().makeShallow();
toc(1)

tic(2)
    sys.setForcingMask();
toc(2)

    if
    (
        SType::nCsComponents == 1
     || sum(sys.diagonal()) == MeshType::numberOfDirections
    )
    {
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
tic(3)
        this->setSingularityConstraint(sys, 0);
toc(3)
        this->solve
        (
            sys,
            this->relTol_,
            this->tolerance_,
            this->minIter_,
            this->maxIter_,
            this->nSweeps_
        );
    }
}

}

}

}
