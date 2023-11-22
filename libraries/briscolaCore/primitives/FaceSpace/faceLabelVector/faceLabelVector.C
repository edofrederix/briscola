#include "faceLabelVector.H"

namespace Foam
{

template<>
const char* const faceLabelVector::csType::typeName = "faceLabelVector";

template<>
const char* const faceLabelVector::csType::componentNames[] =
    {"l", "r", "b", "t", "a", "f"};

template<>
const faceLabelVector faceLabelVector::csType::zero
(
    faceLabelVector::uniform(labelVector::uniform(0))
);

template<>
const faceLabelVector faceLabelVector::csType::one
(
    faceLabelVector::uniform(labelVector::uniform(1))
);

template<>
const faceLabelVector faceLabelVector::csType::max
(
    faceLabelVector::uniform(labelVector::uniform(INT_MAX))
);

template<>
const faceLabelVector faceLabelVector::csType::min
(
    faceLabelVector::uniform(labelVector::uniform(-INT_MAX))
);

template<>
const faceLabelVector faceLabelVector::csType::rootMax
(
    faceLabelVector::uniform(labelVector::uniform(::sqrt(INT_MAX)))
);

template<>
const faceLabelVector faceLabelVector::csType::rootMin
(
    faceLabelVector::uniform(labelVector::uniform(-::sqrt(INT_MAX)))
);

}
