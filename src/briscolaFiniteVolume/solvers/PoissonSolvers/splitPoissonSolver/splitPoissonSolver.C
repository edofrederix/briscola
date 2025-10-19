#include "splitPoissonSolver.H"
#include "exSchemesDivergence.H"
#include "exSchemesFaceGradient.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
splitPoissonSolver<SType,Type,MeshType>::splitPoissonSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,Type,MeshType>(dict,fvMsh),
    solverPtr_
    (
        PoissonSolver<SType,Type,MeshType>::New
        (
            dict.subDict("solver"),
            fvMsh
        )
    ),
    nCorr_(dict.lookupOrDefault<label>("nCorr", 1)),
    prediction_(dict.lookupOrDefault<Switch>("prediction", true))
{
    if (dict.found("preconditioner"))
    {
        preSolverPtr_.reset
        (
            PoissonSolver<SType,Type,MeshType>::New
            (
                dict.subDict("preconditioner"),
                fvMsh
            ).ptr()
        );

        preSolverPtr_->enableFlux();
    }

    if (nCorr_ < 1)
        FatalErrorInFunction
            << "nCorr must be at least 1"
            << endl << abort(FatalError);
}

template<class SType, class Type, class MeshType>
void splitPoissonSolver<SType,Type,MeshType>::solve
(
    meshField<Type,MeshType>& x,
    const meshField<Type,MeshType>* bPtr,
    const faceField<scalar,MeshType>* lambdaPtr,
    const bool ddt,
    const scalar dtFrac
)
{
    if (!lambdaPtr)
    {
        FatalErrorInFunction
            << "Poisson equation has no coefficient so it cannot be split. "
            << endl << abort(FatalError);
    }

    // Set prediction

    if (prediction_)
    {
        const scalar deltaT = x.fvMsh().time().deltaTValue()*dtFrac;
        const scalar deltaT0 = x.fvMsh().time().deltaT0Value();

        x =
            (1.0 + deltaT/deltaT0)*x.oldTime()
          - deltaT/deltaT0*x.oldTime().oldTime();
    }

    // Pre-conditioner

    if (preSolverPtr_.valid())
    {
        preSolverPtr_->solve(x, bPtr, lambdaPtr, ddt, dtFrac);
    }

    // Modified coefficient

    const faceField<scalar,MeshType>& lambda = *lambdaPtr;

    scalar lambda0 = 0.0;

    forAll(lambda, s)
        lambda0 = Foam::max(lambda0, max(gMax(lambda[s])));

    // Modified source

    if (bHatPtr_.empty())
    {
        const word bName
        (
            bPtr ? bPtr->name() : word(x.name() + "_b")
        );

        bHatPtr_.reset
        (
            meshField<Type,MeshType>::New
            (
                bName + "Hat",
                x.fvMsh()
            ).ptr()
        );

        bHatPtr_() = Zero;
    }

    meshField<Type,MeshType>& bHat = bHatPtr_();

    const faceField<scalar,MeshType>& fa =
        x.fvMsh().template metrics<MeshType>().faceAreas();

    for (int corr = 0; corr < nCorr_; corr++)
    {
        if (bPtr)
        {
            bHat = (*bPtr)/lambda0;
        }
        else
        {
            bHat = Zero;
        }

        const faceField<Type,MeshType> phi
        (
            (lambda/lambda0 - 1.0)*ex::faceGrad(x)*fa
        );

        bHat += ex::div(phi);

        solverPtr_->solve(x, &bHat, nullptr, ddt, dtFrac);

        // Compute the flux if needed

        if (this->computeFlux() && corr == nCorr_-1)
        {
            this->initFlux();

            this->fluxPtr_() = solverPtr_->flux() + phi;
            this->fluxPtr_() *= lambda0;
        }
    }
}

}

}

}
