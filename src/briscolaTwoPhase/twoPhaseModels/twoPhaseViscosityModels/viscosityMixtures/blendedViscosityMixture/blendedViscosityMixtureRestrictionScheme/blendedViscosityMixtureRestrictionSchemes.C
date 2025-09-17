#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "blendedViscosityMixtureRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeFaceRestrictionSchemeMeshType(blendedViscosityMixture,faceScalar,colocated);
makeFaceRestrictionSchemeMeshType(blendedViscosityMixture,faceScalar,staggered);

}

}

}
