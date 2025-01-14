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
inline tmp<meshField<FaceSpace<Type>,MeshType>> interp
(
    const meshField<Type,MeshType>& field
)
{
    return interpolationScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        "interp("+field.name()+")"
    )->interp(field);
}

template<class Type, class MeshType>
inline tmp<meshField<FaceSpace<Type>,MeshType>> interp
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<FaceSpace<Type>,MeshType>> tInterp
    (
        interp(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

// Interpolation from colocated to staggered

using ::Foam::briscola::fv::stagInterp;
using ::Foam::briscola::fv::stagFaceInterp;

// Interpolation from staggered to colocated

using ::Foam::briscola::fv::coloInterp;
using ::Foam::briscola::fv::coloFaceInterp;

}

}

}

}

#endif
