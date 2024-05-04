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
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    return faceGradientScheme<Type,MeshType>::New
    (
        "faceGrad("+field.name()+")",
        field.fvMsh()
    )->faceGrad(field);
}

template<class Type, class MeshType>
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> faceGrad
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<meshField<LowerFaceSpace<Type>,MeshType>> tFaceGrad
    (
        faceGrad(tField())
    );

    tField.clear();

    return tFaceGrad;
}

}

}

}

}

#endif
