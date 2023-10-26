#include "geometricVof.H"

#include "addToRunTimeSelectionTable.H"
#include "rectilinearMesh.H"
#include "truncatedHex.H"
#include "truncatedPiped.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(geometricVof, 0);

geometricVof::geometricVof(const dictionary& dict, const fvMesh& fvMsh)
:
    vof(dict, fvMsh),
    rectilinear_(fvMsh.msh()[0].rectilinear() == unitXYZ),
    lve_(rectilinear_),
    normalSchemePtr_
    (
        normalScheme::New
        (
            *this,
            dict.subDict("normalScheme")
        ).ptr()
    )
{}

geometricVof::geometricVof(const geometricVof& vf)
:
    vof(vf),
    rectilinear_(vf.rectilinear_),
    lve_(rectilinear_),
    normalSchemePtr_
    (
        normalScheme::New
        (
            *this,
            dict_.subDict("normalScheme")
        ).ptr()
    )
{}

geometricVof::~geometricVof()
{}

}

}

}
