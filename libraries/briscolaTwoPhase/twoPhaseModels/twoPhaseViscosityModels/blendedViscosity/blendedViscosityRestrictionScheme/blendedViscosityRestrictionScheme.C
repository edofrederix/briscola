#include "blendedViscosityRestrictionScheme.H"

#include "colocated.H"
#include "twoPhaseModel.H"
#include "blendedViscosity.H"
#include "incompressibleTwoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
blendedViscosityRestrictionScheme<MeshType>::blendedViscosityRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<lowerFaceScalar,MeshType>(fvMsh, is)
{
    const twoPhaseModel* modelPtr =
        &fvMsh.db().lookupObject<twoPhaseModel>("briscolaTwoPhaseDict");

    typedef blendedViscosity<incompressibleTwoPhaseModel> modelType;

    if(dynamic_cast<const modelType*>(modelPtr))
    {
        C_ = dynamic_cast<const modelType*>(modelPtr)->C_;
        mu1_ = dynamic_cast<const modelType*>(modelPtr)->mu1_;
        mu2_ = dynamic_cast<const modelType*>(modelPtr)->mu2_;
    }
    else
    {
        FatalErrorInFunction
            << "Two-phase model is of incompatible type"
            << endl << abort(FatalError);
    }
}

}

}

}
