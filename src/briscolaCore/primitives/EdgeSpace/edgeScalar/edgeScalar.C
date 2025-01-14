#include "edgeScalar.H"

namespace Foam
{

template<>
const char* const edgeScalar::csType::typeName = "edgeScalar";

template<>
const char* const edgeScalar::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const labelVector edgeScalar::csType::componentOffsets[] =
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
const edgeScalar edgeScalar::csType::zero(edgeScalar::uniform(0));

template<>
const edgeScalar edgeScalar::csType::one(edgeScalar::uniform(1));

template<>
const edgeScalar edgeScalar::csType::max(edgeScalar::uniform(vGreat));

template<>
const edgeScalar edgeScalar::csType::min(edgeScalar::uniform(-vGreat));

template<>
const edgeScalar edgeScalar::csType::rootMax(edgeScalar::uniform(rootVGreat));

template<>
const edgeScalar edgeScalar::csType::rootMin(edgeScalar::uniform(-rootVGreat));

}