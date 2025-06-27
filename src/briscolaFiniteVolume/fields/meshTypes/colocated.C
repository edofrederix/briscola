#include "colocated.H"
#include "utils.H"

template<>
const int
Foam::briscola::fv::colocated::meshType::numberOfDirections = 1;

template<>
const char* const Foam::briscola::fv::colocated::meshType::directionNames[] =
    {""};

template<>
const Foam::labelVector
Foam::briscola::fv::colocated::meshType::padding[] = {zeroXYZ};

template<>
const Foam::vector
Foam::briscola::fv::colocated::meshType::shift[] = {Foam::vector(0, 0, 0)};

template<>
const char* const
Foam::briscola::fv::colocated::meshType::typeName = "colocated";
