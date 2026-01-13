#include "stencil.H"

namespace Foam
{

template<>
const char* const stencil::csType::typeName = "stencil";

template<>
const char* const stencil::csType::componentNames[] =
    {"c", "l", "r", "b", "t", "a", "f"};

template<>
const labelVector stencil::csType::offsets[] =
    {
        briscola::zeroXYZ,
        briscola::faceOffset0,
        briscola::faceOffset1,
        briscola::faceOffset2,
        briscola::faceOffset3,
        briscola::faceOffset4,
        briscola::faceOffset5
    };

template<>
const stencil stencil::csType::zero(stencil::uniform(0));

template<>
const stencil stencil::csType::one(stencil::uniform(1));

template<>
const stencil stencil::csType::max(stencil::uniform(vGreat));

template<>
const stencil stencil::csType::min(stencil::uniform(-vGreat));

template<>
const stencil stencil::csType::rootMax(stencil::uniform(rootVGreat));

template<>
const stencil stencil::csType::rootMin(stencil::uniform(-rootVGreat));

}