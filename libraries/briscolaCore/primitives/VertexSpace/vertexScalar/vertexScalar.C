#include "vertexScalar.H"

namespace Foam
{

template<>
const char* const vertexScalar::csType::typeName = "vertexScalar";

template<>
const char* const vertexScalar::csType::componentNames[] =
    {"v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7"};

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