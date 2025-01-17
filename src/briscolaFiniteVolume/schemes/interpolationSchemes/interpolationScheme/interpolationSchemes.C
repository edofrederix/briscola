#include "interpolationSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointInterpolationScheme.H"
#include "linearInterpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeInterpolationScheme(scalar,colocated)
makeInterpolationScheme(scalar,staggered)

makeInterpolationScheme(vector,colocated)
makeInterpolationScheme(vector,staggered)

makeInterpolationSchemeType(midPoint,scalar,colocated)
makeInterpolationSchemeType(midPoint,scalar,staggered)

makeInterpolationSchemeType(midPoint,vector,colocated)
makeInterpolationSchemeType(midPoint,vector,staggered)
makeInterpolationSchemeType(linear,scalar,colocated)
makeInterpolationSchemeType(linear,scalar,staggered)

makeInterpolationSchemeType(linear,vector,colocated)
makeInterpolationSchemeType(linear,vector,staggered)

}

}

}
