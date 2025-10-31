#include "defaultPoissonSolver.H"
#include "imSchemes.H"
#include "RungeKuttaScheme.H"
#include "exSchemesFaceGradient.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
defaultPoissonSolver<SType,Type,MeshType>::defaultPoissonSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<SType,Type,MeshType>(dict,fvMsh),
    solverPtr_
    (
        word(dict.lookup("type")) == "default"
      ? solver<SType,Type,MeshType>::New
        (
            dict.subDict("solver"),
            fvMsh
        )
      : solver<SType,Type,MeshType>::New
        (
            dict,
            fvMsh
        )
    )
{}

template<class SType, class Type, class MeshType>
void defaultPoissonSolver<SType,Type,MeshType>::solve
(
    meshField<Type,MeshType>& x,
    const meshField<Type,MeshType>* bPtr,
    const faceField<scalar,MeshType>* lambdaPtr,
    const bool ddt,
    const scalar dtFrac
)
{
    const fvMesh& fvMsh = this->fvMsh_;

    if (sysPtr_.empty() || &sysPtr_->x() != &x)
    {
        sysPtr_.reset
        (
            new linearSystem<SType,Type,MeshType>
            (
                lambdaPtr
              ? word("Poisson("+lambdaPtr->name()+","+x.name()+")")
              : word("Poisson("+x.name()+")"),
                x
            )
        );
    }

    // Initialize the Runge-Kutta stage solution list if we haven't done so yet
    // and if there's a Runge-Kutta scheme pointer

    const RungeKuttaScheme* rkSchemePtr = nullptr;

    if (fvMsh.db().template foundObject<RungeKuttaScheme>("rkScheme"))
        rkSchemePtr =
            &fvMsh.db().template lookupObject<RungeKuttaScheme>("rkScheme");

    if (rkSchemePtr && !this->stages_.size())
    {
        this->stages_.resize(rkSchemePtr->nStages());

        forAll(this->stages_, i)
            this->stages_.set
            (
                i,
                meshField<Type,MeshType>::New
                (
                    x.name() + "_" + Foam::name(i),
                    fvMsh
                )
            );
    }

    linearSystem<SType,Type,MeshType>& sys = sysPtr_();

    sys = im::laplacian<SType>(lambdaPtr,x);

    if (bPtr)
        sys += (*bPtr);

    if (ddt)
        sys -= im::ddt(x)/dtFrac;

    const bool constMatrix = !ddt && !lambdaPtr;

    // If we have a Runge-Kutta scheme pointer, use the previous stage solution
    // as initial guess and store the solution for the next iteration. This
    // improves convergence.

    if (rkSchemePtr)
        x = this->stages_[rkSchemePtr->stage()-1];

    solverPtr_->solve(sys,constMatrix);

    if (rkSchemePtr)
        this->stages_[rkSchemePtr->stage()-1] = x;

    // Compute the flux if needed

    if (this->computeFlux())
    {
        this->initFlux();

        const faceField<scalar,MeshType>& fa =
            x.fvMsh().template metrics<MeshType>().faceAreas();

        if (lambdaPtr)
        {
            this->fluxPtr_() = (*lambdaPtr)*ex::faceGrad(x)*fa;

            this->fluxPtr_->rename
            (
                lambdaPtr->name()
              + "*faceGrad(" + x.name() + ")"
              + "*" + fa.name()
            );
        }
        else
        {
            this->fluxPtr_() = ex::faceGrad(x)*fa;

            this->fluxPtr_->rename
            (
                "*faceGrad(" + x.name() + ")"
              + "*" + fa.name()
            );
        }
    }
}

}

}

}
