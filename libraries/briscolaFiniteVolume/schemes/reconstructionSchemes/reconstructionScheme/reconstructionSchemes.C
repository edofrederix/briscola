#include "reconstructionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointReconstructionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeReconstructionScheme(scalar)
makeReconstructionScheme(vector)

makeReconstructionSchemeType(midPoint,scalar)
makeReconstructionSchemeType(midPoint,vector)

}

}

}
