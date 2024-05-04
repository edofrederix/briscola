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

template<class Type, class MeshType>
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<meshField<LowerFaceSpace<Type>,MeshType>> tInterp
    (
        interp(tField())
    );

    tField.clear();

    return tInterp;
}

// Interpolation from colocated to staggered

using ::Foam::briscola::fv::stagInterp;
using ::Foam::briscola::fv::stagFaceInterp;

}

}

}

}

#endif
