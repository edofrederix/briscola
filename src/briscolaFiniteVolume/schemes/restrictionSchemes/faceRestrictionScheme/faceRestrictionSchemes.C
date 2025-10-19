#include "faceRestrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "areaWeightedFaceRestrictionScheme.H"
#include "averageFaceRestrictionScheme.H"
#include "fluxFaceRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Face restriction schemes

makeFaceRestrictionSchemeBase(scalar,colocated,"areaWeighted");
makeFaceRestrictionSchemeBase(vector,colocated,"areaWeighted");

makeFaceRestrictionSchemeBase(scalar,staggered,"areaWeighted");
makeFaceRestrictionSchemeBase(vector,staggered,"areaWeighted");

makeFaceRestrictionScheme(average,scalar,colocated);
makeFaceRestrictionScheme(average,vector,colocated);
makeFaceRestrictionScheme(average,scalar,staggered);
makeFaceRestrictionScheme(average,vector,staggered);

makeFaceRestrictionScheme(areaWeighted,scalar,colocated);
makeFaceRestrictionScheme(areaWeighted,vector,colocated);
makeFaceRestrictionScheme(areaWeighted,scalar,staggered);
makeFaceRestrictionScheme(areaWeighted,vector,staggered);

makeFaceRestrictionScheme(flux,scalar,colocated);
makeFaceRestrictionScheme(flux,vector,colocated);

}

}

}
