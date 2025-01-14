#include "splitPoissonSolver.H"
#include "exSchemes.H"

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
    )
{}

template<class SType, class Type, class MeshType>
void splitPoissonSolver<SType,Type,MeshType>::solve
(
    meshField<Type,MeshType>& x,
    const meshField<Type,MeshType>* bPtr,
    const meshField<faceScalar,MeshType>* lambdaPtr,
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

    // Modified coefficient

    const meshField<faceScalar,MeshType>& lambda = *lambdaPtr;

    const faceScalar lambda0f = max(gMax(lambda));
    const scalar lambda0 =
        Foam::max
        (
            Foam::max
            (
                lambda0f.left(),
                lambda0f.bottom()
            ),
            lambda0f.aft()
        );

    // Modified source

    if (bHatPtr_.empty())
    {
        const word bName
        (
            bPtr ? bPtr->name() : word(x.name() + "_b")
        );

        bHatPtr_.reset
        (
            new meshField<Type,MeshType>
            (
                bName + "Hat",
                x.fvMsh()
            )
        );

        bHatPtr_() = Zero;
    }

    meshField<Type,MeshType>& bHat = bHatPtr_();

    if (bPtr)
    {
        bHat = (*bPtr);
        bHat /= lambda0;
    }
    else
    {
        bHat = Zero;
    }

    const scalar deltaT = x.fvMsh().time().deltaTValue()*dtFrac;
    const scalar deltaT0 = x.fvMsh().time().deltaT0Value();

    const meshField<Type,MeshType> xHat
    (
        (1.0 + deltaT/deltaT0)*x.oldTime()
      - deltaT/deltaT0*x.oldTime().oldTime()
    );

    bHat += ex::laplacian(lambda/lambda0 - 1.0, xHat);

    // Solve approximate constant coefficient system following Dodd & Ferrante
    // (2014), Eq. (16)

    solverPtr_->solve(x, &bHat, nullptr, ddt, dtFrac);

    // Compute the flux if needed

    if (this->computeFlux())
    {
        this->initFlux();

        const meshField<faceScalar,MeshType>& fa =
            x.fvMsh().template metrics<MeshType>().faceAreas();

        // Dodd & Ferrante (2014), Eq. (14)

        this->fluxPtr_() =
            lambda0*solverPtr_->flux()
          + (lambda - lambda0)*ex::faceGrad(xHat)*fa;

        this->fluxPtr_->rename
        (
            lambda.name() + "*" + solverPtr_->flux().name()
        );
    }
}

}

}

}
