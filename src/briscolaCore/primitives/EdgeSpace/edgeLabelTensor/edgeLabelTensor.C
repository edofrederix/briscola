#include "edgeLabelTensor.H"

namespace Foam
{

template<>
const char* const edgeLabelTensor::csType::typeName = "edgeLabelTensor";

template<>
const char* const edgeLabelTensor::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const labelVector edgeLabelTensor::csType::offsets[] =
    {
        briscola::edgeOffset0,
        briscola::edgeOffset1,
        briscola::edgeOffset2,
        briscola::edgeOffset3,
        briscola::edgeOffset4,
        briscola::edgeOffset5,
        briscola::edgeOffset6,
        briscola::edgeOffset7,
        briscola::edgeOffset8,
        briscola::edgeOffset9,
        briscola::edgeOffset10,
        briscola::edgeOffset11
    };

template<>
const edgeLabelTensor edgeLabelTensor::csType::zero
(
    edgeLabelTensor::uniform(labelTensor::uniform(0))
);

template<>
const edgeLabelTensor edgeLabelTensor::csType::one
(
    edgeLabelTensor::uniform(labelTensor::uniform(1))
);

template<>
const edgeLabelTensor edgeLabelTensor::csType::max
(
    edgeLabelTensor::uniform(labelTensor::uniform(INT_MAX))
);

template<>
const edgeLabelTensor edgeLabelTensor::csType::min
(
    edgeLabelTensor::uniform(labelTensor::uniform(-INT_MAX))
);

template<>
const edgeLabelTensor edgeLabelTensor::csType::rootMax
(
    edgeLabelTensor::uniform(labelTensor::uniform(::sqrt(INT_MAX)))
);

template<>
const edgeLabelTensor edgeLabelTensor::csType::rootMin
(
    edgeLabelTensor::uniform(labelTensor::uniform(-::sqrt(INT_MAX)))
);

}
