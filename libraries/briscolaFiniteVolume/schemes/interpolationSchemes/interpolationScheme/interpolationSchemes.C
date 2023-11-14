#include "interpolationSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointInterpolationScheme.H"
#include "curvatureInterpolationScheme.H"

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

makeInterpolationSchemeType(curvature,scalar,colocated)
makeInterpolationSchemeType(curvature,scalar,staggered)

makeInterpolationSchemeType(curvature,vector,colocated)
makeInterpolationSchemeType(curvature,vector,staggered)

}

}

}
