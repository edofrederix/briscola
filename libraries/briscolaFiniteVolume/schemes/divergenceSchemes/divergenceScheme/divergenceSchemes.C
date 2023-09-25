#include "divergenceSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointGaussDivergenceScheme.H"
#include "upwindDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeDivergenceScheme(scalar,colocated)
makeDivergenceScheme(scalar,staggered)

makeDivergenceScheme(vector,colocated)
makeDivergenceScheme(vector,staggered)

makeDivergenceSchemeType(midPointGauss,scalar,colocated)
makeDivergenceSchemeType(midPointGauss,scalar,staggered)

makeDivergenceSchemeType(midPointGauss,vector,colocated)
makeDivergenceSchemeType(midPointGauss,vector,staggered)

makeDivergenceSchemeType(upwind,scalar,colocated)
makeDivergenceSchemeType(upwind,scalar,staggered)

makeDivergenceSchemeType(upwind,vector,colocated)
makeDivergenceSchemeType(upwind,vector,staggered)

}

}

}
