#include "edgeVector.H"

namespace Foam
{

template<>
const char* const edgeVector::csType::typeName = "edgeVector";

template<>
const char* const edgeVector::csType::componentNames[] =
    {"ba", "ta", "bf", "tf", "la", "ra", "lf", "rf", "lb", "rb", "lt", "rt"};

template<>
const labelVector edgeVector::csType::offsets[] =
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
const edgeVector edgeVector::csType::zero
(
    edgeVector::uniform(vector::uniform(0))
);

template<>
const edgeVector edgeVector::csType::one
(
    edgeVector::uniform(vector::uniform(1))
);

template<>
const edgeVector edgeVector::csType::max
(
    edgeVector::uniform(vector::uniform(vGreat))
);

template<>
const edgeVector edgeVector::csType::min
(
    edgeVector::uniform(vector::uniform(-vGreat))
);

template<>
const edgeVector edgeVector::csType::rootMax
(
    edgeVector::uniform(vector::uniform(rootVGreat))
);

template<>
const edgeVector edgeVector::csType::rootMin
(
    edgeVector::uniform(vector::uniform(-rootVGreat))
);

}