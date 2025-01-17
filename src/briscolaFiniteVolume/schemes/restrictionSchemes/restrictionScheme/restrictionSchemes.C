#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "averageRestrictionScheme.H"
#include "linearRestrictionScheme.H"
#include "volumeWeightedRestrictionScheme.H"
#include "harmonicVolumeWeightedRestrictionScheme.H"

#include "coloFaceAverageRestrictionScheme.H"
#include "coloFaceAreaWeightedRestrictionScheme.H"
#include "coloHarmonicFaceAreaWeightedRestrictionScheme.H"
#include "coloFluxRestrictionScheme.H"

#include "stagFaceAverageRestrictionScheme.H"

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

makeRestrictionSchemeType
(
    coloFaceAverage,
    faceAverage,
    faceScalar,
    colocated
);
makeRestrictionSchemeType
(
    coloFaceAverage,
    faceAverage,
    faceVector,
    colocated
);

makeRestrictionSchemeType
(
    stagFaceAverage,
    faceAverage,
    faceScalar,
    staggered
);
makeRestrictionSchemeType
(
    stagFaceAverage,
    faceAverage,
    faceVector,
    staggered
);

makeRestrictionSchemeType
(
    coloFaceAreaWeighted,
    faceAreaWeighted,
    faceScalar,
    colocated
);
makeRestrictionSchemeType
(
    coloFaceAreaWeighted,
    faceAreaWeighted,
    faceVector,
    colocated
);

makeRestrictionSchemeType
(
    coloHarmonicFaceAreaWeighted,
    harmonicFaceAreaWeighted,
    faceScalar,
    colocated
);

makeRestrictionSchemeNoTemplate(coloFlux,faceScalar,colocated);

}

}

}
