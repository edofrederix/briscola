#include "faceVector.H"

namespace Foam
{

template<>
const char* const faceVector::csType::typeName = "faceVector";

template<>
const char* const faceVector::csType::componentNames[] =
    {"l", "r", "b", "t", "a", "f"};

template<>
const labelVector faceVector::csType::offsets[] =
    {
        briscola::faceOffset0,
        briscola::faceOffset1,
        briscola::faceOffset2,
        briscola::faceOffset3,
        briscola::faceOffset4,
        briscola::faceOffset5
    };

template<>
const faceVector faceVector::csType::zero
(
    faceVector::uniform(vector::uniform(0))
);

template<>
const faceVector faceVector::csType::one
(
    faceVector::uniform(vector::uniform(1))
);

template<>
const faceVector faceVector::csType::max
(
    faceVector::uniform(vector::uniform(vGreat))
);

template<>
const faceVector faceVector::csType::min
(
    faceVector::uniform(vector::uniform(-vGreat))
);

template<>
const faceVector faceVector::csType::rootMax
(
    faceVector::uniform(vector::uniform(rootVGreat))
);

template<>
const faceVector faceVector::csType::rootMin
(
    faceVector::uniform(vector::uniform(-rootVGreat))
);

}