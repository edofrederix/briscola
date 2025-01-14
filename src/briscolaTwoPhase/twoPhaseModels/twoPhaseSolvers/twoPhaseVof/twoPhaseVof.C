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
    vfPtr_(vof::New(*this, dict.subDict("vof")))
{}

template<class BaseModel>
twoPhaseVof<BaseModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    BaseModel(tpm),
    vfPtr_(tpm.vf().clone())
{
    this->normalSchemePtr_->correct();
    BaseModel::correctMixture();
}

template<class BaseModel>
twoPhaseVof<BaseModel>::~twoPhaseVof()
{}

template<class BaseModel>
void twoPhaseVof<BaseModel>::correct()
{
    vfPtr_->solve(this->coloFaceFlux()());
    BaseModel::correctMixture();
}

}

}

}
