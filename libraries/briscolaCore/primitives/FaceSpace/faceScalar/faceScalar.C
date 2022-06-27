#include "faceScalar.H"

namespace Foam
{

template<>
const char* const faceScalar::csType::typeName = "faceScalar";

template<>
const char* const faceScalar::csType::componentNames[] = {"l", "r", "b", "t", "a", "f"};

template<>
const faceScalar faceScalar::csType::zero(faceScalar::uniform(0));

template<>
const faceScalar faceScalar::csType::one(faceScalar::uniform(1));

template<>
const faceScalar faceScalar::csType::max(faceScalar::uniform(vGreat));

template<>
const faceScalar faceScalar::csType::min(faceScalar::uniform(-vGreat));

template<>
const faceScalar faceScalar::csType::rootMax(faceScalar::uniform(rootVGreat));

template<>
const faceScalar faceScalar::csType::rootMin(faceScalar::uniform(-rootVGreat));

}