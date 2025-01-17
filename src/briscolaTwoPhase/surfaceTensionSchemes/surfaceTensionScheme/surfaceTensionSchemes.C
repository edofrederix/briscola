#include "surfaceTensionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "constantSigma.H"
#include "zeroSigma.H"
#include "Brackbill.H"
#include "noSurfaceTension.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// No surface tension

makeSurfaceTensionScheme
(
    none,
    noSurfaceTension,
    zeroSigma
)

// Brackbill with constant sigma value

makeSurfaceTensionScheme
(
    Brackbill,
    Brackbill,
    constantSigma
)

}

}

}
