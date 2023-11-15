#include "staggered.H"
#include "briscolaUtils.H"

template<>
const int
Foam::briscola::fv::staggered::meshType::numberOfDirections = 3;

template<>
const char* const Foam::briscola::fv::staggered::meshType::directionNames[] =
    {"x", "y", "z"};

template<>
const Foam::labelVector
Foam::briscola::fv::staggered::meshType::padding[] = {unitX, unitY, unitZ};

template<>
const Foam::vector
Foam::briscola::fv::staggered::meshType::shift[] =
    {
        Foam::vector(-0.5, 0, 0),
        Foam::vector(0, -0.5, 0),
        Foam::vector(0, 0, -0.5)
    };

template<>
const char* const
Foam::briscola::fv::staggered::meshType::typeName = "staggered";