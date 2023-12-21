#ifndef exSchemesGradient_H
#define exSchemesGradient_H

#include "gradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Centered gradient

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<Type,MeshType>& field
)
{
    return gradientScheme<Type,MeshType>::New
    (
        "grad("+field.name()+")",
        field.fvMsh()
    )->grad(field);
}

// Centered gradient from a face field

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<FaceSpace<Type>,MeshType>& field
)
{
    return centerGrad(field);
}

// Staggered gradient from a colocated field

template<class Type>
inline tmp<meshField<Type,staggered>>
stagGrad
(
    const meshField<Type,colocated>& field
)
{
    return gradientScheme<Type,colocated>::New
    (
        "stagGrad("+field.name()+")",
        field.fvMsh()
    )->stagGrad(field);
}


}

}

}

}

#endif
