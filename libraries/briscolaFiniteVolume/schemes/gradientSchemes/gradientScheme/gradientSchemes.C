#include "gradientSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointGaussGradientScheme.H"

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

}

}

}
