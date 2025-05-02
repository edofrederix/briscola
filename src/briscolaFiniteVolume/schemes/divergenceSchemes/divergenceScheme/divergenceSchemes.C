#include "divergenceSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointGaussDivergenceScheme.H"
#include "linearGaussDivergenceScheme.H"
#include "limitedGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeDivergenceScheme(stencil,scalar,colocated)
makeDivergenceScheme(stencil,scalar,staggered)

makeDivergenceScheme(stencil,vector,colocated)
makeDivergenceScheme(stencil,vector,staggered)


makeDivergenceSchemeType(midPointGauss,stencil,scalar,colocated)
makeDivergenceSchemeType(midPointGauss,stencil,scalar,staggered)

makeDivergenceSchemeType(midPointGauss,stencil,vector,colocated)
makeDivergenceSchemeType(midPointGauss,stencil,vector,staggered)


makeDivergenceSchemeType(linearGauss,stencil,scalar,colocated)
makeDivergenceSchemeType(linearGauss,stencil,scalar,staggered)

makeDivergenceSchemeType(linearGauss,stencil,vector,colocated)
makeDivergenceSchemeType(linearGauss,stencil,vector,staggered)


makeDivergenceSchemeType(limitedGauss,stencil,scalar,colocated)
makeDivergenceSchemeType(limitedGauss,stencil,scalar,staggered)

makeDivergenceSchemeType(limitedGauss,stencil,vector,colocated)
makeDivergenceSchemeType(limitedGauss,stencil,vector,staggered)

}

}

}
