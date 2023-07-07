#include "nearestPointInterpolator.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(nearestPointInterpolator, 0);

addToRunTimeSelectionTable
(
    pointInterpolator,
    nearestPointInterpolator,
    dictionary
);

nearestPointInterpolator::nearestPointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points
)
:
    pointInterpolator(fvMsh, points)
{}

nearestPointInterpolator::nearestPointInterpolator
(
    const nearestPointInterpolator& interp
)
:
    pointInterpolator(interp)
{}

#define INTERPFUNC(TYPE)                                                    \
List<TYPE> nearestPointInterpolator::operator()                             \
(                                                                           \
    const meshField<TYPE,colocated>& field                                  \
)                                                                           \
{                                                                           \
    return this->interp(field);                                             \
}

INTERPFUNC(scalar);
INTERPFUNC(vector);
INTERPFUNC(tensor);
INTERPFUNC(sphericalTensor);
INTERPFUNC(symmTensor);
INTERPFUNC(diagTensor);
INTERPFUNC(faceScalar);
INTERPFUNC(faceVector);

#undef INTERPFUNC

template<class Type>
List<Type> nearestPointInterpolator::interp
(
    const meshField<Type,colocated>& field
)
{
    const meshDirection<Type,colocated>& f = field[0][0];

    List<Type> values(points_.size());

    forAll(points_, i)
    {
        if (indices_[i] != -unitXYZ)
        {
            values[i] = f(indices_[i]);
        }
        else
        {
            values[i] = Zero;
        }
    }

    this->gatherScatter(values);

    return values;
}

}

}

}
