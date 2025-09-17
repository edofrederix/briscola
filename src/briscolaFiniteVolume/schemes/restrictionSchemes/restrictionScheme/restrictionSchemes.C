#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "averageRestrictionScheme.H"
#include "linearRestrictionScheme.H"
#include "volumeWeightedRestrictionScheme.H"
#include "harmonicVolumeWeightedRestrictionScheme.H"

#include "harmonicFaceAreaWeightedRestrictionScheme.H"
#include "faceAreaWeightedRestrictionScheme.H"
#include "faceAverageRestrictionScheme.H"
#include "fluxRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Cell restriction schemes

makeRestrictionSchemeBase(scalar,colocated,"linear");
makeRestrictionSchemeBase(scalar,staggered,"linear");
makeRestrictionSchemeBase(vector,colocated,"linear");
makeRestrictionSchemeBase(vector,staggered,"linear");
makeRestrictionSchemeBase(tensor,colocated,"linear");
makeRestrictionSchemeBase(tensor,staggered,"linear");

makeRestrictionScheme(average,scalar,colocated);
makeRestrictionScheme(average,scalar,staggered);
makeRestrictionScheme(average,vector,colocated);
makeRestrictionScheme(average,vector,staggered);
makeRestrictionScheme(average,tensor,colocated);
makeRestrictionScheme(average,tensor,staggered);

makeRestrictionScheme(linear,scalar,colocated);
makeRestrictionScheme(linear,scalar,staggered);
makeRestrictionScheme(linear,vector,colocated);
makeRestrictionScheme(linear,vector,staggered);
makeRestrictionScheme(linear,tensor,colocated);
makeRestrictionScheme(linear,tensor,staggered);

makeRestrictionScheme(volumeWeighted,scalar,colocated);
makeRestrictionScheme(volumeWeighted,scalar,staggered);
makeRestrictionScheme(volumeWeighted,vector,colocated);
makeRestrictionScheme(volumeWeighted,vector,staggered);
makeRestrictionScheme(volumeWeighted,tensor,colocated);
makeRestrictionScheme(volumeWeighted,tensor,staggered);

makeRestrictionScheme(harmonicVolumeWeighted,scalar,colocated);
makeRestrictionScheme(harmonicVolumeWeighted,scalar,staggered);

// Face restriction schemes

makeRestrictionSchemeBase(faceScalar,colocated,"faceAverage");
makeRestrictionSchemeBase(faceVector,colocated,"faceAverage");

makeRestrictionSchemeBase(faceScalar,staggered,"faceAverage");
makeRestrictionSchemeBase(faceVector,staggered,"faceAverage");

makeFaceRestrictionScheme(faceAverage,scalar,faceScalar,colocated);
makeFaceRestrictionScheme(faceAverage,vector,faceVector,colocated);
makeFaceRestrictionScheme(faceAverage,scalar,faceScalar,staggered);
makeFaceRestrictionScheme(faceAverage,vector,faceVector,staggered);

makeFaceRestrictionScheme(faceAreaWeighted,scalar,faceScalar,colocated);
makeFaceRestrictionScheme(faceAreaWeighted,vector,faceVector,colocated);
makeFaceRestrictionScheme(faceAreaWeighted,scalar,faceScalar,staggered);
makeFaceRestrictionScheme(faceAreaWeighted,vector,faceVector,staggered);

makeFaceRestrictionScheme(harmonicFaceAreaWeighted,scalar,faceScalar,colocated);
makeFaceRestrictionScheme(harmonicFaceAreaWeighted,scalar,faceScalar,staggered);

makeFaceRestrictionScheme(flux,scalar,faceScalar,colocated);
makeFaceRestrictionScheme(flux,vector,faceVector,colocated);

}

}

}
