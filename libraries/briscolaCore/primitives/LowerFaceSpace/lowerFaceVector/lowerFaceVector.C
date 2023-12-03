#include "lowerFaceVector.H"

namespace Foam
{

template<>
const char* const lowerFaceVector::csType::typeName = "lowerFaceVector";

template<>
const char* const lowerFaceVector::csType::componentNames[] =
    {"l", "b", "a"};

template<>
const labelVector lowerFaceVector::csType::componentOffsets[] =
    {
        briscola::faceOffset0,
        briscola::faceOffset2,
        briscola::faceOffset4
    };

template<>
const lowerFaceVector lowerFaceVector::csType::zero
(
    lowerFaceVector::uniform(vector::uniform(0))
);

template<>
const lowerFaceVector lowerFaceVector::csType::one
(
    lowerFaceVector::uniform(vector::uniform(1))
);

template<>
const lowerFaceVector lowerFaceVector::csType::max
(
    lowerFaceVector::uniform(vector::uniform(vGreat))
);

template<>
const lowerFaceVector lowerFaceVector::csType::min
(
    lowerFaceVector::uniform(vector::uniform(-vGreat))
);

template<>
const lowerFaceVector lowerFaceVector::csType::rootMax
(
    lowerFaceVector::uniform(vector::uniform(rootVGreat))
);

template<>
const lowerFaceVector lowerFaceVector::csType::rootMin
(
    lowerFaceVector::uniform(vector::uniform(-rootVGreat))
);

}