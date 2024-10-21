#include "blendedViscosityRestrictionScheme.H"

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
{}

}

}

}
