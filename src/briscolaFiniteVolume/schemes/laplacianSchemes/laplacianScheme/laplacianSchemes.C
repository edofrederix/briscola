#include "laplacianSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearGaussLaplacianScheme.H"
#include "stencilLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLaplacianScheme(stencil,scalar,colocated)
makeLaplacianScheme(stencil,vector,colocated)
makeLaplacianScheme(stencil,scalar,staggered)
makeLaplacianScheme(stencil,vector,staggered)

makeLaplacianSchemeType(linearGauss,stencil,scalar,colocated)
makeLaplacianSchemeType(linearGauss,stencil,vector,colocated)
makeLaplacianSchemeType(linearGauss,stencil,scalar,staggered)
makeLaplacianSchemeType(linearGauss,stencil,vector,staggered)

addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,scalar,colocated)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,vector,colocated)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,scalar,staggered)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,vector,staggered)

}

}

}
