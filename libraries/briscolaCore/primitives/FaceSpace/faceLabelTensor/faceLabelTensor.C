#include "faceLabelTensor.H"

namespace Foam
{

template<>
const char* const faceLabelTensor::csType::typeName = "faceLabelTensor";

template<>
const char* const faceLabelTensor::csType::componentNames[] =
    {"l", "r", "b", "t", "a", "f"};

template<>
const labelVector faceLabelTensor::csType::componentOffsets[] =
    {
        briscola::faceOffset0,
        briscola::faceOffset1,
        briscola::faceOffset2,
        briscola::faceOffset3,
        briscola::faceOffset4,
        briscola::faceOffset5
    };

template<>
const faceLabelTensor faceLabelTensor::csType::zero
(
    faceLabelTensor::uniform(labelTensor::uniform(0))
);

template<>
const faceLabelTensor faceLabelTensor::csType::one
(
    faceLabelTensor::uniform(labelTensor::uniform(1))
);

template<>
const faceLabelTensor faceLabelTensor::csType::max
(
    faceLabelTensor::uniform(labelTensor::uniform(INT_MAX))
);

template<>
const faceLabelTensor faceLabelTensor::csType::min
(
    faceLabelTensor::uniform(labelTensor::uniform(-INT_MAX))
);

template<>
const faceLabelTensor faceLabelTensor::csType::rootMax
(
    faceLabelTensor::uniform(labelTensor::uniform(::sqrt(INT_MAX)))
);

template<>
const faceLabelTensor faceLabelTensor::csType::rootMin
(
    faceLabelTensor::uniform(labelTensor::uniform(-::sqrt(INT_MAX)))
);

}
