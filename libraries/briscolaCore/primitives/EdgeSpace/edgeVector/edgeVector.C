#include "edgeVector.H"

namespace Foam
{

template<>
const char* const edgeVector::csType::typeName = "edgeVector";

template<>
const char* const edgeVector::csType::componentNames[] =
    {"x0", "x1", "x2", "x3", "y0", "y1", "y2", "y3", "z0", "z1", "z2", "z3"};

template<>
const edgeVector edgeVector::csType::zero(edgeVector::uniform(vector::uniform(0)));

template<>
const edgeVector edgeVector::csType::one(edgeVector::uniform(vector::uniform(1)));

template<>
const edgeVector edgeVector::csType::max(edgeVector::uniform(vector::uniform(vGreat)));

template<>
const edgeVector edgeVector::csType::min(edgeVector::uniform(vector::uniform(-vGreat)));

template<>
const edgeVector edgeVector::csType::rootMax(edgeVector::uniform(vector::uniform(rootVGreat)));

template<>
const edgeVector edgeVector::csType::rootMin(edgeVector::uniform(vector::uniform(-rootVGreat)));

}