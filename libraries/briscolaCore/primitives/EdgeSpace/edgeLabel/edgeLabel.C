#include "edgeLabel.H"

namespace Foam
{

template<>
const char* const edgeLabel::csType::typeName = "edgeLabel";

template<>
const char* const edgeLabel::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const edgeLabel edgeLabel::csType::zero(edgeLabel::uniform(0));

template<>
const edgeLabel edgeLabel::csType::one(edgeLabel::uniform(1));

template<>
const edgeLabel edgeLabel::csType::max(edgeLabel::uniform(pTraits<label>::max));

template<>
const edgeLabel edgeLabel::csType::min(edgeLabel::uniform(pTraits<label>::min));

template<>
const edgeLabel edgeLabel::csType::rootMax
(
    edgeLabel::uniform(pTraits<label>::rootMax)
);

template<>
const edgeLabel edgeLabel::csType::rootMin
(
    edgeLabel::uniform(pTraits<label>::rootMin)
);

}