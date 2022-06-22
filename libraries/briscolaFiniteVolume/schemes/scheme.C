#include "scheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(scheme, 0);

scheme::scheme(const dictionary& dict, const fvMesh& fvMsh)
:
    dict_(dict),
    fvMsh_(fvMsh)
{}

scheme::scheme(const scheme& s)
:
    dict_(s.dict_),
    fvMsh_(s.fvMsh_)
{}

scheme::~scheme()
{}

}

}

}
