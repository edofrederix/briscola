#include "edgeLabelVector.H"

namespace Foam
{

template<>
const char* const edgeLabelVector::csType::typeName = "edgeLabelVector";

template<>
const char* const edgeLabelVector::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const labelVector edgeLabelVector::csType::componentOffsets[] =
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
