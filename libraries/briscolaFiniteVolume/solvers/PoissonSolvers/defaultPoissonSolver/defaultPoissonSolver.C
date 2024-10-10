#include "defaultPoissonSolver.H"
#include "imSchemes.H"

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
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    const bool ddt,
    const scalar dtFrac
)
{
    if (sysPtr_.empty() || &sysPtr_->x() != &x)
    {
        sysPtr_.reset
        (
            new linearSystem<SType,Type,MeshType>(x)
        );
    }

    linearSystem<SType,Type,MeshType>& sys = sysPtr_();

    sys = im::laplacian<SType>(lambdaPtr,x);

    if (bPtr)
        sys += (*bPtr);

    if (ddt)
        sys -= im::ddt(x, 1.0/dtFrac);

    const bool constMatrix = !ddt && !lambdaPtr;

    solverPtr_->solve(sys,constMatrix);

    // Compute the flux if needed

    if (this->computeFlux())
    {
        this->initFlux();

        const meshField<faceScalar,MeshType>& fa =
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
