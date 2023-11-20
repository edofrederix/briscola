#include "vertexLabel.H"

namespace Foam
{

template<>
const char* const vertexLabel::csType::typeName = "vertexLabel";

template<>
const char* const vertexLabel::csType::componentNames[] =
    {"lba", "rba", "lta", "rta", "lbf", "rbf", "ltf", "rtf"};

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
    vertexLabel::uniform(pTraits<label>::rootMax)
);

template<>
const vertexLabel vertexLabel::csType::rootMin
(
    vertexLabel::uniform(pTraits<label>::rootMin)
);

}