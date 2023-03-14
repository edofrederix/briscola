#include "vertexVector.H"

namespace Foam
{

template<>
const char* const vertexVector::csType::typeName = "vertexVector";

template<>
const char* const vertexVector::csType::componentNames[] =
    {"v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7"};

template<>
const vertexVector vertexVector::csType::zero(vertexVector::uniform(vector::uniform(0)));

template<>
const vertexVector vertexVector::csType::one(vertexVector::uniform(vector::uniform(1)));

template<>
const vertexVector vertexVector::csType::max(vertexVector::uniform(vector::uniform(vGreat)));

template<>
const vertexVector vertexVector::csType::min(vertexVector::uniform(vector::uniform(-vGreat)));

template<>
const vertexVector vertexVector::csType::rootMax(vertexVector::uniform(vector::uniform(rootVGreat)));

template<>
const vertexVector vertexVector::csType::rootMin(vertexVector::uniform(vector::uniform(-rootVGreat)));

}