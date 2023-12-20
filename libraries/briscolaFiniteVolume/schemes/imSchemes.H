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

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>> div
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    return divergenceScheme<SType,Type,MeshType>::New
    (
        "div("+phi.name()+","+field.name()+")",
        field.fvMsh()
    )->imDiv(phi,field);
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>> laplacian
(
    const meshField<lowerFaceScalar,MeshType>& lambda,
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<SType,Type,MeshType>::New
    (
        "laplacian("+lambda.name()+","+field.name()+")",
        field.fvMsh()
    )->imLaplacian(lambda,field);
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>> laplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<SType,Type,MeshType>::New
    (
        "laplacian("
      + (lambdaPtr ? lambdaPtr->name()+"," : word(""))
      + field.name()
      + ")",
        field.fvMsh()
    )->imLaplacian(lambdaPtr,field);
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>> laplacian
(
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<SType,Type,MeshType>::New
    (
        "laplacian("+field.name()+")",
        field.fvMsh()
    )->imLaplacian(field);
}

template<class Type, class MeshType>
inline tmp<linearSystem<diagStencil,Type,MeshType>> source
(
    const meshField<scalar,MeshType>& coeff,
    meshField<Type,MeshType>& field
)
{
    const_cast<meshField<scalar,MeshType>&>(coeff).restrict();

    tmp<linearSystem<diagStencil,Type,MeshType>> tSys
    (
        new linearSystem<diagStencil,Type,MeshType>(field)
    );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    meshField<diagStencil,MeshType>& A = Sys.A();

    A = Zero;

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) = cv(l,d,i,j,k)*coeff(l,d,i,j,k);

    Sys.b() = Zero;

    const_cast<meshField<scalar,MeshType>&>(coeff).makeShallow();

    return tSys;
}

// Default stencil types

// Full stencil for any divergence operator

template<class Type, class MeshType>
inline tmp<linearSystem<stencil,Type,MeshType>> div
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    return div<stencil>(phi,field);
}

// Symmetric stencil for Laplacians on colocated meshes

template<class Type>
inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
(
    const meshField<lowerFaceScalar,colocated>& lambda,
    meshField<Type,colocated>& field
)
{
    return laplacian<symmStencil>(lambda,field);
}

template<class Type>
inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
(
    const meshField<lowerFaceScalar,colocated>* lambdaPtr,
    meshField<Type,colocated>& field
)
{
    return laplacian<symmStencil>(lambdaPtr,field);
}

template<class Type>
inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
(
    meshField<Type,colocated>& field
)
{
    return laplacian<symmStencil>(field);
}

// Full stencil for Laplacians on staggered meshes

template<class Type>
inline tmp<linearSystem<stencil,Type,staggered>> laplacian
(
    const meshField<lowerFaceScalar,staggered>& lambda,
    meshField<Type,staggered>& field
)
{
    return laplacian<stencil>(lambda,field);
}

template<class Type>
inline tmp<linearSystem<stencil,Type,staggered>> laplacian
(
    const meshField<lowerFaceScalar,staggered>* lambdaPtr,
    meshField<Type,staggered>& field
)
{
    return laplacian<stencil>(lambdaPtr,field);
}

template<class Type>
inline tmp<linearSystem<stencil,Type,staggered>> laplacian
(
    meshField<Type,staggered>& field
)
{
    return laplacian<stencil>(field);
}

}

}

}

}

#endif
