#include "laplacianSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearGaussLaplacianScheme.H"
#include "coloLinearGaussLaplacianScheme.H"
#include "stagLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLaplacianScheme(symmStencil,scalar,colocated)
makeLaplacianScheme(symmStencil,vector,colocated)
makeLaplacianScheme(stencil,scalar,staggered)
makeLaplacianScheme(stencil,vector,staggered)

makeLaplacianSchemeType(linearGauss,symmStencil,scalar,colocated)
makeLaplacianSchemeType(linearGauss,symmStencil,vector,colocated)
makeLaplacianSchemeType(linearGauss,stencil,scalar,staggered)
makeLaplacianSchemeType(linearGauss,stencil,vector,staggered)

addSpecificLaplacianSchemeType(coloLinearGauss,symmStencil,scalar,colocated)
addSpecificLaplacianSchemeType(coloLinearGauss,symmStencil,vector,colocated)
addSpecificLaplacianSchemeType(stagLinearGauss,stencil,scalar,staggered)
addSpecificLaplacianSchemeType(stagLinearGauss,stencil,vector,staggered)

}

}

}
