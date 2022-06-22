#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "averageRestrictionScheme.H"
#include "linearRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionScheme(scalar,colocated);
makeRestrictionScheme(scalar,staggered);
makeRestrictionScheme(vector,colocated);
makeRestrictionScheme(vector,staggered);
makeRestrictionScheme(tensor,colocated);
makeRestrictionScheme(tensor,staggered);

makeRestrictionSchemeType(average,scalar,colocated);
makeRestrictionSchemeType(average,scalar,staggered);
makeRestrictionSchemeType(average,vector,colocated);
makeRestrictionSchemeType(average,vector,staggered);
makeRestrictionSchemeType(average,tensor,colocated);
makeRestrictionSchemeType(average,tensor,staggered);

makeRestrictionSchemeType(linear,scalar,colocated);
makeRestrictionSchemeType(linear,scalar,staggered);
makeRestrictionSchemeType(linear,vector,colocated);
makeRestrictionSchemeType(linear,vector,staggered);
makeRestrictionSchemeType(linear,tensor,colocated);
makeRestrictionSchemeType(linear,tensor,staggered);

}

}

}
