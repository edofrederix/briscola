#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "averageRestrictionScheme.H"
#include "linearRestrictionScheme.H"
#include "fluxRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionScheme(scalar,colocated,"linear");
makeRestrictionScheme(scalar,staggered,"linear");
makeRestrictionScheme(vector,colocated,"linear");
makeRestrictionScheme(vector,staggered,"linear");
makeRestrictionScheme(tensor,colocated,"linear");
makeRestrictionScheme(tensor,staggered,"linear");
makeRestrictionScheme(faceScalar,colocated,"linear");
makeRestrictionScheme(faceScalar,staggered,"linear");

makeRestrictionSchemeType(average,scalar,colocated);
makeRestrictionSchemeType(average,scalar,staggered);
makeRestrictionSchemeType(average,vector,colocated);
makeRestrictionSchemeType(average,vector,staggered);
makeRestrictionSchemeType(average,tensor,colocated);
makeRestrictionSchemeType(average,tensor,staggered);
makeRestrictionSchemeType(average,faceScalar,colocated);
makeRestrictionSchemeType(average,faceScalar,staggered);

makeRestrictionSchemeType(linear,scalar,colocated);
makeRestrictionSchemeType(linear,scalar,staggered);
makeRestrictionSchemeType(linear,vector,colocated);
makeRestrictionSchemeType(linear,vector,staggered);
makeRestrictionSchemeType(linear,tensor,colocated);
makeRestrictionSchemeType(linear,tensor,staggered);
makeRestrictionSchemeType(linear,faceScalar,colocated);
makeRestrictionSchemeType(linear,faceScalar,staggered);

makeRestrictionSchemeSingleType(flux,faceScalar,colocated);
makeRestrictionSchemeSingleType(flux,faceScalar,staggered);

}

}

}
