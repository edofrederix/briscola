#include "symmStencil.H"

namespace Foam
{

template<>
const char* const symmStencil::csType::typeName = "symmStencil";

template<>
const char* const symmStencil::csType::componentNames[] =
    {"c", "l", "b", "a"};

template<>
const labelVector symmStencil::csType::componentOffsets[] =
    {
        briscola::zeroXYZ,
        briscola::faceOffset0,
        briscola::faceOffset2,
        briscola::faceOffset4
    };

template<>
const symmStencil symmStencil::csType::zero(symmStencil::uniform(0));

template<>
const symmStencil symmStencil::csType::one(symmStencil::uniform(1));

template<>
const symmStencil symmStencil::csType::max(symmStencil::uniform(vGreat));

template<>
const symmStencil symmStencil::csType::min(symmStencil::uniform(-vGreat));

template<>
const symmStencil symmStencil::csType::rootMax
(
    symmStencil::uniform(rootVGreat)
);

template<>
const symmStencil symmStencil::csType::rootMin
(
    symmStencil::uniform(-rootVGreat)
);

}