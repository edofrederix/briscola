#include "solverFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<>
const char* NamedEnum<smootherType,3>::names[] =
{
    "RBGS",
    "LEXGS",
    "JAC"
};

}

}

}
