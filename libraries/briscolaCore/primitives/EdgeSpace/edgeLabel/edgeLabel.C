#include "edgeLabel.H"

namespace Foam
{

template<>
const char* const edgeLabel::csType::typeName = "edgeLabel";

template<>
const char* const edgeLabel::csType::componentNames[] =
    {"x0", "x1", "x2", "x3", "y0", "y1", "y2", "y3", "z0", "z1", "z2", "z3"};

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