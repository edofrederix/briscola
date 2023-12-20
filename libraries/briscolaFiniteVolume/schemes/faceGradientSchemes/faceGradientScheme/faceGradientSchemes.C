#include "faceGradientSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearFaceGradientScheme.H"

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

makeFaceGradientSchemeType(linear,scalar,colocated)
makeFaceGradientSchemeType(linear,scalar,staggered)

makeFaceGradientSchemeType(linear,vector,colocated)
makeFaceGradientSchemeType(linear,vector,staggered)

}

}

}
