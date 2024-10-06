#include "faceFluxSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointFaceFluxScheme.H"
#include "linearFaceFluxScheme.H"

// Also compile non-templated faceFluxScheme
#include "faceFluxScheme.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(faceFluxScheme, 0);
defineRunTimeSelectionTable(faceFluxScheme, Istream);

makeFaceFluxSchemeType(midPoint)
makeFaceFluxSchemeType(linear)

}

}

}
