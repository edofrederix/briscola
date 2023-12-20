#include "lowerFaceScalar.H"

namespace Foam
{

template<>
const char* const lowerFaceScalar::csType::typeName = "lowerFaceScalar";

template<>
const char* const lowerFaceScalar::csType::componentNames[] =
    {"l", "b", "a"};

template<>
const labelVector lowerFaceScalar::csType::componentOffsets[] =
    {
        briscola::faceOffset0,
        briscola::faceOffset2,
        briscola::faceOffset4
    };

template<>
const lowerFaceScalar lowerFaceScalar::csType::zero
(
    lowerFaceScalar::uniform(0)
);

template<>
const lowerFaceScalar lowerFaceScalar::csType::one
(
    lowerFaceScalar::uniform(1)
);

template<>
const lowerFaceScalar lowerFaceScalar::csType::max
(
    lowerFaceScalar::uniform(vGreat)
);

template<>
const lowerFaceScalar lowerFaceScalar::csType::min
(
    lowerFaceScalar::uniform(-vGreat)
);

template<>
const lowerFaceScalar lowerFaceScalar::csType::rootMax
(
    lowerFaceScalar::uniform(rootVGreat)
);

template<>
const lowerFaceScalar lowerFaceScalar::csType::rootMin
(
    lowerFaceScalar::uniform(-rootVGreat)
);

}
