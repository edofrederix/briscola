#include "vertexLabelTensor.H"

namespace Foam
{

template<>
const char* const vertexLabelTensor::csType::typeName = "vertexLabelTensor";

template<>
const char* const vertexLabelTensor::csType::componentNames[] =
    {"lba", "rba", "lta", "rta", "lbf", "rbf", "ltf", "rtf"};

template<>
const vertexLabelTensor vertexLabelTensor::csType::zero
(
    vertexLabelTensor::uniform(labelTensor::uniform(0))
);

template<>
const vertexLabelTensor vertexLabelTensor::csType::one
(
    vertexLabelTensor::uniform(labelTensor::uniform(1))
);

template<>
const vertexLabelTensor vertexLabelTensor::csType::max
(
    vertexLabelTensor::uniform(labelTensor::uniform(INT_MAX))
);

template<>
const vertexLabelTensor vertexLabelTensor::csType::min
(
    vertexLabelTensor::uniform(labelTensor::uniform(-INT_MAX))
);

template<>
const vertexLabelTensor vertexLabelTensor::csType::rootMax
(
    vertexLabelTensor::uniform(labelTensor::uniform(::sqrt(INT_MAX)))
);

template<>
const vertexLabelTensor vertexLabelTensor::csType::rootMin
(
    vertexLabelTensor::uniform(labelTensor::uniform(-::sqrt(INT_MAX)))
);

}
