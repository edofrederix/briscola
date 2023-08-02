#include "vof.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vof, 0);

vof::vof(const IOdictionary& dict, const fvMesh& fvMsh)
:
    regIOobject(dict, true),
    fvMsh_(fvMsh),
    alpha_
    (
        "alpha",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    ),
    lve_(*this),
    normalSchemePtr_
    (
        normalScheme::New
        (
            *this,
            dict.subDict("normalScheme")
        ).ptr()
    )
{}

vof::~vof()
{}

}

}

}
