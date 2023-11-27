#include "faceLabel.H"

namespace Foam
{

template<>
const char* const faceLabel::csType::typeName = "faceLabel";

template<>
const char* const faceLabel::csType::componentNames[] =
    {"l", "r", "b", "t", "a", "f"};

template<>
const labelVector faceLabel::csType::componentOffsets[] =
    {
        briscola::faceOffset0,
        briscola::faceOffset1,
        briscola::faceOffset2,
        briscola::faceOffset3,
        briscola::faceOffset4,
        briscola::faceOffset5
    };

template<>
const faceLabel faceLabel::csType::zero(faceLabel::uniform(0));

template<>
const faceLabel faceLabel::csType::one(faceLabel::uniform(1));

template<>
const faceLabel faceLabel::csType::max(faceLabel::uniform(INT_MAX));

template<>
const faceLabel faceLabel::csType::min(faceLabel::uniform(-INT_MAX));

template<>
const faceLabel faceLabel::csType::rootMax(faceLabel::uniform(::sqrt(INT_MAX)));

template<>
const faceLabel faceLabel::csType::rootMin
(
    faceLabel::uniform(-::sqrt(INT_MAX))
);

}