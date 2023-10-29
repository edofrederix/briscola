#include "surfaceTensionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "constantSigma.H"
#include "Brackbill.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

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
