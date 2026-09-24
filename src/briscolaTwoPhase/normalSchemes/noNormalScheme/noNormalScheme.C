#include "noNormalScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(noNormalScheme, 0);
addToRunTimeSelectionTable(normalScheme, noNormalScheme, dictionary);

noNormalScheme::noNormalScheme
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dict, alpha)
{}

noNormalScheme::noNormalScheme
(
    const fvMesh& fvMsh,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dictionary::null, alpha)
{}

noNormalScheme::noNormalScheme(const noNormalScheme& s)
:
    normalScheme(s)
{}

noNormalScheme::~noNormalScheme()
{}

void noNormalScheme::correct()
{}

}

}

}
