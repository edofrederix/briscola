#include "faceGradientSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointFaceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeFaceGradientScheme(scalar,colocated)
makeFaceGradientScheme(scalar,staggered)

makeFaceGradientScheme(vector,colocated)
makeFaceGradientScheme(vector,staggered)

makeFaceGradientSchemeType(midPoint,scalar,colocated)
makeFaceGradientSchemeType(midPoint,scalar,staggered)

makeFaceGradientSchemeType(midPoint,vector,colocated)
makeFaceGradientSchemeType(midPoint,vector,staggered)

}

}

}
