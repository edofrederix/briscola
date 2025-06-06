#include "vertexLabel.H"

namespace Foam
{

template<>
const char* const vertexLabel::csType::typeName = "vertexLabel";

template<>
const char* const vertexLabel::csType::componentNames[] =
    {"lba", "rba", "lta", "rta", "lbf", "rbf", "ltf", "rtf"};

template<>
const labelVector vertexLabel::csType::componentOffsets[] =
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
const vertexLabel vertexLabel::csType::zero(vertexLabel::uniform(0));

template<>
const vertexLabel vertexLabel::csType::one(vertexLabel::uniform(1));

template<>
const vertexLabel vertexLabel::csType::max
(
    vertexLabel::uniform(pTraits<label>::max)
);

template<>
const vertexLabel vertexLabel::csType::min
(
    vertexLabel::uniform(pTraits<label>::min)
);

template<>
const vertexLabel vertexLabel::csType::rootMax
(
    vertexLabel::uniform(::sqrt(INT_MAX))
);

template<>
const vertexLabel vertexLabel::csType::rootMin
(
    vertexLabel::uniform(-::sqrt(INT_MAX))
);

}