#include "scheme.H"

#include "faceFields.H"
#include "meshFields.H"
#include "meshLevels.H"
#include "meshDirections.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(scheme, 0);

dummyIstream scheme::nullStream = dummyIstream();

scheme::scheme(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh)
{}

scheme::scheme(const scheme& s)
:
    fvMsh_(s.fvMsh_)
{}

scheme::~scheme()
{}

Istream& scheme::getStream
(
    const fvMesh& fvMsh,
    const word& dictName,
    const word& schemeName
)
{
    const dictionary& dict = fvMsh.schemeDict().subDict(dictName);

    if (dict.found(schemeName))
    {
        return dict.lookup(schemeName);
    }
    else if (dict.found("default"))
    {
        return dict.lookup("default");
    }
    else
    {
        FatalErrorInFunction
            << "Scheme " << schemeName << " or default scheme not found in "
            << "dictionary " << dictName
            << endl << abort(FatalError);

        return dict.lookup("none");
    }
}

}

}

}
