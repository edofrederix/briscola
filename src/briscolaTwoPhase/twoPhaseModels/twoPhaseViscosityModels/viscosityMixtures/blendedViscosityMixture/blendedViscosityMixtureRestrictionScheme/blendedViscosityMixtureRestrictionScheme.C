#include "blendedViscosityMixtureRestrictionScheme.H"
#include "twoPhaseModel.H"

#ifdef NoRepository
#undef NoRepository
#endif

#include "blendedViscosityMixture.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
blendedViscosityMixtureRestrictionScheme<MeshType>::
blendedViscosityMixtureRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<faceScalar,MeshType>(fvMsh, is)
{
    const twoPhaseModel& model =
        fvMsh.db().lookupObject<twoPhaseModel>("briscolaTwoPhaseDict");

    this->C_ = model.dict().lookupOrDefault<scalar>("C", 8);

    // This currently assumes a Newtonian fluid
    mu1_ = readScalar(model.dict().lookup("mu1"));
    mu2_ = readScalar(model.dict().lookup("mu2"));
}

}

}

}
