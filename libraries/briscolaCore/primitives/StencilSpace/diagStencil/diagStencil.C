#include "diagStencil.H"

namespace Foam
{

template<>
const char* const diagStencil::csType::typeName = "diagStencil";

template<>
const char* const diagStencil::csType::componentNames[] = {"c"};

template<>
const labelVector diagStencil::csType::componentOffsets[] =
    {
        briscola::zeroXYZ
    };

template<>
const diagStencil diagStencil::csType::zero(diagStencil::uniform(0));

template<>
const diagStencil diagStencil::csType::one(diagStencil::uniform(1));

template<>
const diagStencil diagStencil::csType::max(diagStencil::uniform(vGreat));

template<>
const diagStencil diagStencil::csType::min(diagStencil::uniform(-vGreat));

template<>
const diagStencil diagStencil::csType::rootMax(diagStencil::uniform(rootVGreat));

template<>
const diagStencil diagStencil::csType::rootMin(diagStencil::uniform(-rootVGreat));

}
