#include "nearestPointInterpolator.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebugWithName
(
    nearestPointInterpolator<colocated>,
    "nearest",
    0
);

defineTemplateTypeNameAndDebugWithName
(
    nearestPointInterpolator<staggered>,
    "nearest",
    0
);

addTemplatedToRunTimeSelectionTable
(
    pointInterpolator,
    nearestPointInterpolator,
    colocated,
    dictionary
);

addTemplatedToRunTimeSelectionTable
(
    pointInterpolator,
    nearestPointInterpolator,
    staggered,
    dictionary
);

template<class MeshType>
nearestPointInterpolator<MeshType>::nearestPointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points,
    const bool global,
    const label l,
    const label d
)
:
    pointInterpolator<MeshType>(fvMsh, points, global, l, d)
{}

template<class MeshType>
nearestPointInterpolator<MeshType>::nearestPointInterpolator
(
    const nearestPointInterpolator<MeshType>& interp
)
:
    pointInterpolator<MeshType>(interp)
{}

template<class MeshType>
nearestPointInterpolator<MeshType>::~nearestPointInterpolator()
{}

#define INTERPFUNC(TYPE)                                                    \
template<class MeshType>                                                    \
List<TYPE> nearestPointInterpolator<MeshType>::operator()                   \
(                                                                           \
    const meshField<TYPE,MeshType>& field                                   \
) const                                                                     \
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

template<class MeshType>
template<class Type>
List<Type> nearestPointInterpolator<MeshType>::interp
(
    const meshField<Type,MeshType>& field
) const
{
    List<Type> values(this->points_.size(), Zero);

    const label me = Pstream::myProcNo();

    forAll(this->points_, i)
        if (this->providers_[i] == me)
            values[i] = field(this->l_, this->d_, this->indices_[i]);

    return this->combine(values);
}

// Instantiate

template class nearestPointInterpolator<colocated>;
template class nearestPointInterpolator<staggered>;

}

}

}
