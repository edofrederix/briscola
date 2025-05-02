#include "limiterSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "twoPhaseLimiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeLimiterSchemeType(twoPhase,scalar,colocated);
makeLimiterSchemeType(twoPhase,scalar,staggered);
makeLimiterSchemeType(twoPhase,vector,colocated);
makeLimiterSchemeType(twoPhase,vector,staggered);

}

}

}
