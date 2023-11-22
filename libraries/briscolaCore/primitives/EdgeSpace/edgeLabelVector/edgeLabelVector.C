#include "edgeLabelVector.H"

namespace Foam
{

template<>
const char* const edgeLabelVector::csType::typeName = "edgeLabelVector";

template<>
const char* const edgeLabelVector::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const edgeLabelVector edgeLabelVector::csType::zero
(
    edgeLabelVector::uniform(labelVector::uniform(0))
);

template<>
const edgeLabelVector edgeLabelVector::csType::one
(
    edgeLabelVector::uniform(labelVector::uniform(1))
);

template<>
const edgeLabelVector edgeLabelVector::csType::max
(
    edgeLabelVector::uniform(labelVector::uniform(INT_MAX))
);

template<>
const edgeLabelVector edgeLabelVector::csType::min
(
    edgeLabelVector::uniform(labelVector::uniform(-INT_MAX))
);

template<>
const edgeLabelVector edgeLabelVector::csType::rootMax
(
    edgeLabelVector::uniform(labelVector::uniform(::sqrt(INT_MAX)))
);

template<>
const edgeLabelVector edgeLabelVector::csType::rootMin
(
    edgeLabelVector::uniform(labelVector::uniform(-::sqrt(INT_MAX)))
);

}
