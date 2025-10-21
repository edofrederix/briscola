#include "edgeLabel.H"

namespace Foam
{

template<>
const char* const edgeLabel::csType::typeName = "edgeLabel";

template<>
const char* const edgeLabel::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const labelVector edgeLabel::csType::offsets[] =
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
const edgeLabel edgeLabel::csType::zero(edgeLabel::uniform(0));

template<>
const edgeLabel edgeLabel::csType::one(edgeLabel::uniform(1));

template<>
const edgeLabel edgeLabel::csType::max(edgeLabel::uniform(pTraits<label>::max));

template<>
const edgeLabel edgeLabel::csType::min(edgeLabel::uniform(pTraits<label>::min));

template<>
const edgeLabel edgeLabel::csType::rootMax(edgeLabel::uniform(::sqrt(INT_MAX)));

template<>
const edgeLabel edgeLabel::csType::rootMin
(
    edgeLabel::uniform(-::sqrt(INT_MAX))
);

}