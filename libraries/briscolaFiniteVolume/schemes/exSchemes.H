#ifndef exSchemes_H
#define exSchemes_H

#include "ddtScheme.H"
#include "divergenceScheme.H"
#include "faceFluxScheme.H"
#include "faceGradientScheme.H"
#include "gradientScheme.H"
#include "interpolationScheme.H"
#include "laplacianScheme.H"
#include "reconstructionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> div
(
    const meshField<Hex<Type>,MeshType>& phi
)
{
    return explicitDiv(phi);
}

template<class Type>
inline tmp<meshField<Type,colocated>> coloDiv
(
    const meshField<Type,staggered>& field
)
{
    return explicitColoDiv(field);
}

inline tmp<meshField<hexScalar,colocated>> faceFlux
(
    const meshField<vector,colocated>& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

inline tmp<meshField<hexScalar,staggered>> faceFlux
(
    const meshField<scalar,staggered>& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

template<class Type, class MeshType>
inline tmp<meshField<Hex<Type>,MeshType>> faceGrad
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

template<class Type, class MeshType>
inline tmp<meshField<Hex<Type>,MeshType>> interp
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

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<Type,staggered>& field
)
{
    return reconstructionScheme<Type>::New
    (
        "reconstruct("+field.name()+")",
        field.fvMsh()
    )->reconstruct(field);
}

}

}

}

}

#endif
