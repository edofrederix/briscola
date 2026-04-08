#ifndef exSchemesFaceGradient_H
#define exSchemesFaceGradient_H

#include "faceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Face gradient

template<class Type, class MeshType>
inline tmp<faceField<Type,MeshType>> faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    return faceGradientScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        "faceGrad("+field.name()+")"
    )->faceGrad(field);
}

template<class Type, class MeshType>
inline tmp<faceField<Type,MeshType>> faceGrad
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    if (tField.isTmp())
        const_cast<tmp<meshField<Type,MeshType>>&>(tField)
            ->correctBoundaryConditions();

    tmp<faceField<Type,MeshType>> tFaceGrad =
        faceGrad(tField());

    if (tField.isTmp())
        tField.clear();

    return tFaceGrad;
}

}

}

}

}

#endif
