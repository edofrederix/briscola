#ifndef imSchemes_H
#define imSchemes_H

#include "ddtScheme.H"
#include "divergenceScheme.H"
#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace im
{

template<class Type, class MeshType>
inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
(
    meshField<Type,MeshType>& field
)
{
    return ddtScheme<Type,MeshType>::New
    (
        "ddt("+field.name()+")",
        field.fvMsh()
    )->ddt(field);
}

template<class Type, class MeshType>
inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
(
    const meshField<scalar,MeshType>& coeff,
    meshField<Type,MeshType>& field
)
{
    return ddtScheme<Type,MeshType>::New
    (
        "ddt("+coeff.name()+","+field.name()+")",
        field.fvMsh()
    )->ddt(coeff,field);
}

template<class Type, class MeshType>
inline tmp<linearSystem<stencil,Type,MeshType>> div
(
    const meshField<hexScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    return divergenceScheme<Type,MeshType>::New
    (
        "div("+phi.name()+","+field.name()+")",
        field.fvMsh()
    )->div(phi,field);
}

template<class Type, class MeshType>
inline tmp<linearSystem<stencil,Type,MeshType>> laplacian
(
    const meshField<hexScalar,MeshType>& lambda,
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<Type,MeshType>::New
    (
        "laplacian("+lambda.name()+","+field.name()+")",
        field.fvMsh()
    )->laplacian(lambda,field);
}

template<class Type, class MeshType>
inline tmp<linearSystem<stencil,Type,MeshType>> laplacian
(
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<Type,MeshType>::New
    (
        "laplacian("+field.name()+")",
        field.fvMsh()
    )->laplacian(field);
}

}

}

}

}

#endif
