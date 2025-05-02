#include "limiterSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "upwindLimiterScheme.H"
#include "vanLeerLimiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLimiterScheme(scalar,colocated);
makeLimiterScheme(scalar,staggered);
makeLimiterScheme(vector,colocated);
makeLimiterScheme(vector,staggered);

makeLimiterSchemeType(upwind,scalar,colocated);
makeLimiterSchemeType(upwind,scalar,staggered);
makeLimiterSchemeType(upwind,vector,colocated);
makeLimiterSchemeType(upwind,vector,staggered);

makeLimiterSchemeType(vanLeer,scalar,colocated);
makeLimiterSchemeType(vanLeer,scalar,staggered);
makeLimiterSchemeType(vanLeer,vector,colocated);
makeLimiterSchemeType(vanLeer,vector,staggered);

}

}

}
