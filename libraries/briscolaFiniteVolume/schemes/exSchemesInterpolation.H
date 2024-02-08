#ifndef exSchemesInterpolation_H
#define exSchemesInterpolation_H

#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Interpolation to the faces

template<class Type, class MeshType>
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
(
    const meshField<Type,MeshType>& field
)
{
    return interpolationScheme<Type,MeshType>::New
    (
        "interp("+field.name()+")",
        field.fvMsh()
    )->interp(field);
}

}

}

}

}

#endif
