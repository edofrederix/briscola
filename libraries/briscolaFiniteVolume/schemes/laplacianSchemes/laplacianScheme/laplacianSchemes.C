#include "laplacianSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLaplacianScheme(scalar,colocated)
makeLaplacianScheme(scalar,staggered)

makeLaplacianScheme(vector,colocated)
makeLaplacianScheme(vector,staggered)

makeLaplacianSchemeType(linearGauss,scalar,colocated)
makeLaplacianSchemeType(linearGauss,scalar,staggered)

makeLaplacianSchemeType(linearGauss,vector,colocated)
makeLaplacianSchemeType(linearGauss,vector,staggered)

}

}

}
