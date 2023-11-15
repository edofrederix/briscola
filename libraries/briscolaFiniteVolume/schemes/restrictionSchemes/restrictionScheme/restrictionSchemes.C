#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "averageRestrictionScheme.H"
#include "linearRestrictionScheme.H"
#include "volumeWeightedRestrictionScheme.H"
#include "harmonicVolumeWeightedRestrictionScheme.H"

#include "faceAverageRestrictionScheme.H"
#include "faceAreaWeightedRestrictionScheme.H"
#include "harmonicFaceAreaWeightedRestrictionScheme.H"
#include "fluxRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Cell restriction schemes

makeRestrictionScheme(scalar,colocated,"linear");
makeRestrictionScheme(scalar,staggered,"linear");
makeRestrictionScheme(vector,colocated,"linear");
makeRestrictionScheme(vector,staggered,"linear");
makeRestrictionScheme(tensor,colocated,"linear");
makeRestrictionScheme(tensor,staggered,"linear");

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

makeRestrictionSchemeType(volumeWeighted,scalar,colocated);
makeRestrictionSchemeType(volumeWeighted,scalar,staggered);
makeRestrictionSchemeType(volumeWeighted,vector,colocated);
makeRestrictionSchemeType(volumeWeighted,vector,staggered);
makeRestrictionSchemeType(volumeWeighted,tensor,colocated);
makeRestrictionSchemeType(volumeWeighted,tensor,staggered);

makeRestrictionSchemeType(harmonicVolumeWeighted,scalar,colocated);
makeRestrictionSchemeType(harmonicVolumeWeighted,scalar,staggered);

// Face restriction schemes

makeRestrictionScheme(faceScalar,colocated,"faceAverage");
makeRestrictionScheme(faceScalar,staggered,"faceAverage");
makeRestrictionScheme(faceVector,colocated,"faceAverage");
makeRestrictionScheme(faceVector,staggered,"faceAverage");

makeRestrictionSchemeType(faceAverage,faceScalar,colocated);
makeRestrictionSchemeType(faceAverage,faceScalar,staggered);
makeRestrictionSchemeType(faceAverage,faceVector,colocated);
makeRestrictionSchemeType(faceAverage,faceVector,staggered);

makeRestrictionSchemeType(faceAreaWeighted,faceScalar,colocated);
makeRestrictionSchemeType(faceAreaWeighted,faceScalar,staggered);
makeRestrictionSchemeType(faceAreaWeighted,faceVector,colocated);
makeRestrictionSchemeType(faceAreaWeighted,faceVector,staggered);

makeRestrictionSchemeType(harmonicFaceAreaWeighted,faceScalar,colocated);
makeRestrictionSchemeType(harmonicFaceAreaWeighted,faceScalar,staggered);

makeRestrictionSchemeSingleType(flux,faceScalar,colocated);
makeRestrictionSchemeSingleType(flux,faceScalar,staggered);

}

}

}
