#include "vertexLabelVector.H"

namespace Foam
{

template<>
const char* const vertexLabelVector::csType::typeName = "vertexLabelVector";

template<>
const char* const vertexLabelVector::csType::componentNames[] =
    {"lba", "rba", "lta", "rta", "lbf", "rbf", "ltf", "rtf"};

template<>
const labelVector vertexLabelVector::csType::componentOffsets[] =
    {
        briscola::vertexOffset0,
        briscola::vertexOffset1,
        briscola::vertexOffset2,
        briscola::vertexOffset3,
        briscola::vertexOffset4,
        briscola::vertexOffset5,
        briscola::vertexOffset6,
        briscola::vertexOffset7,
    };

template<>
const vertexLabelVector vertexLabelVector::csType::zero
(
    vertexLabelVector::uniform(labelVector::uniform(0))
);

template<>
const vertexLabelVector vertexLabelVector::csType::one
(
    vertexLabelVector::uniform(labelVector::uniform(1))
);

template<>
const vertexLabelVector vertexLabelVector::csType::max
(
    vertexLabelVector::uniform(labelVector::uniform(INT_MAX))
);

template<>
const vertexLabelVector vertexLabelVector::csType::min
(
    vertexLabelVector::uniform(labelVector::uniform(-INT_MAX))
);

template<>
const vertexLabelVector vertexLabelVector::csType::rootMax
(
    vertexLabelVector::uniform(labelVector::uniform(::sqrt(INT_MAX)))
);

template<>
const vertexLabelVector vertexLabelVector::csType::rootMin
(
    vertexLabelVector::uniform(labelVector::uniform(-::sqrt(INT_MAX)))
);

}
