#include "faceFluxSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(faceFluxScheme, 0);
defineRunTimeSelectionTable(faceFluxScheme, dictionary);

makeFaceFluxSchemeType(midPoint)

}

}

}
