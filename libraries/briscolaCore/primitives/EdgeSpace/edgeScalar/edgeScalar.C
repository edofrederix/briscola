#include "edgeScalar.H"

namespace Foam
{

template<>
const char* const edgeScalar::csType::typeName = "edgeScalar";

template<>
const char* const edgeScalar::csType::componentNames[] =
    {"x0", "x1", "x2", "x3", "y0", "y1", "y2", "y3", "z0", "z1", "z2", "z3"};

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