#include "gradientSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointGaussGradientScheme.H"
#include "linearGaussGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeGradientScheme(scalar,colocated)
makeGradientScheme(scalar,staggered)

makeGradientScheme(vector,colocated)
makeGradientScheme(vector,staggered)

makeGradientSchemeType(midPointGauss,scalar,colocated)
makeGradientSchemeType(midPointGauss,scalar,staggered)

makeGradientSchemeType(midPointGauss,vector,colocated)
makeGradientSchemeType(midPointGauss,vector,staggered)

makeGradientSchemeType(linearGauss,scalar,colocated)
makeGradientSchemeType(linearGauss,scalar,staggered)

makeGradientSchemeType(linearGauss,vector,colocated)
makeGradientSchemeType(linearGauss,vector,staggered)

}

}

}
