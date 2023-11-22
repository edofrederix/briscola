#include "vertexScalar.H"

namespace Foam
{

template<>
const char* const vertexScalar::csType::typeName = "vertexScalar";

template<>
const char* const vertexScalar::csType::componentNames[] =
    {"lba", "rba", "lta", "rta", "lbf", "rbf", "ltf", "rtf"};

template<>
const vertexScalar vertexScalar::csType::zero(vertexScalar::uniform(0));

template<>
const vertexScalar vertexScalar::csType::one(vertexScalar::uniform(1));

template<>
const vertexScalar vertexScalar::csType::max(vertexScalar::uniform(vGreat));

template<>
const vertexScalar vertexScalar::csType::min(vertexScalar::uniform(-vGreat));

template<>
const vertexScalar vertexScalar::csType::rootMax(vertexScalar::uniform(rootVGreat));

template<>
const vertexScalar vertexScalar::csType::rootMin(vertexScalar::uniform(-rootVGreat));

}