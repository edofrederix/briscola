#include "faceDotGradientSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearFaceDotGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeFaceDotGradientScheme(vector,colocated)
makeFaceDotGradientScheme(scalar,staggered)

makeFaceDotGradientSchemeType(linear,vector,colocated)
makeFaceDotGradientSchemeType(linear,scalar,staggered)

}

}

}
