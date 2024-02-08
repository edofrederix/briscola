#include "laplacianSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearGaussLaplacianScheme.H"
#include "symmStencilLinearGaussLaplacianScheme.H"
#include "stencilLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLaplacianScheme(symmStencil,scalar,colocated)
makeLaplacianScheme(symmStencil,vector,colocated)
makeLaplacianScheme(symmStencil,scalar,staggered)
makeLaplacianScheme(symmStencil,vector,staggered)

makeLaplacianScheme(stencil,scalar,colocated)
makeLaplacianScheme(stencil,vector,colocated)
makeLaplacianScheme(stencil,scalar,staggered)
makeLaplacianScheme(stencil,vector,staggered)

makeLaplacianSchemeType(linearGauss,symmStencil,scalar,colocated)
makeLaplacianSchemeType(linearGauss,symmStencil,vector,colocated)
makeLaplacianSchemeType(linearGauss,symmStencil,scalar,staggered)
makeLaplacianSchemeType(linearGauss,symmStencil,vector,staggered)

makeLaplacianSchemeType(linearGauss,stencil,scalar,colocated)
makeLaplacianSchemeType(linearGauss,stencil,vector,colocated)
makeLaplacianSchemeType(linearGauss,stencil,scalar,staggered)
makeLaplacianSchemeType(linearGauss,stencil,vector,staggered)

addSpecificLaplacianSchemeType(symmStencilLinearGauss,symmStencil,scalar,colocated)
addSpecificLaplacianSchemeType(symmStencilLinearGauss,symmStencil,vector,colocated)
addSpecificLaplacianSchemeType(symmStencilLinearGauss,symmStencil,scalar,staggered)
addSpecificLaplacianSchemeType(symmStencilLinearGauss,symmStencil,vector,staggered)

addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,scalar,colocated)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,vector,colocated)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,scalar,staggered)
addSpecificLaplacianSchemeType(stencilLinearGauss,stencil,vector,staggered)

}

}

}
