#include "twoPhaseVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
twoPhaseVof<BaseModel>::twoPhaseVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    normalSchemePtr_
    (
        normalScheme::New(*this, dict.subDict("normalScheme")).ptr()
    ),
    surfaceTensionSchemePtr_
    (
        surfaceTensionScheme::New
        (
            fvMsh,
            dict.subDict("surfaceTensionScheme"),
            normalSchemePtr_(),
            this->alpha()
        ).ptr()
    ),
    vfPtr_
    (
        vof::New
        (
            fvMsh,
            dict.subDict("vof"),
            normalSchemePtr_(),
            this->alpha()
        )
    )
{}

template<class BaseModel>
twoPhaseVof<BaseModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    BaseModel(tpm),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false),
    vfPtr_(tpm.vf().clone())
{
    normalSchemePtr_->correct();
    surfaceTensionSchemePtr_->correct();

    BaseModel::correctMixture();
}

template<class BaseModel>
twoPhaseVof<BaseModel>::~twoPhaseVof()
{}

template<class BaseModel>
void twoPhaseVof<BaseModel>::correct()
{
    vfPtr_->solve(this->coloFaceFlux()());
    surfaceTensionSchemePtr_->correct();

    BaseModel::correctMixture();
}

}

}

}
