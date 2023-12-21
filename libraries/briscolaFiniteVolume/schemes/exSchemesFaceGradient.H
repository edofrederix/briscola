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

}

}

}

}

#endif
