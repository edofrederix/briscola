#include "boundaryCondition.H"

namespace Foam
{

template<>
const char* NamedEnum<briscola::fv::boundaryConditionBaseType,7>::names[] =
{
    "dummy",
    "empty",
    "parallel",
    "periodic",
    "Dirichlet",
    "Neumann",
    "Robin"
};

}
