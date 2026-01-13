#include "MG.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<>
const char* NamedEnum<MGCycleType,3>::names[] =
{
    "V",
    "W",
    "F"
};

template<>
const char* NamedEnum<MGCoarseMode,2>::names[] =
{
    "smooth",
    "direct"
};

}

}

}
