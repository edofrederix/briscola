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
    const meshField<faceScalar,MeshType>& phi,
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
    const meshField<faceScalar,MeshType>& lambda,
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

template<class Type, class MeshType>
inline tmp<linearSystem<diagStencil,Type,MeshType>> source
(
    const meshField<scalar,MeshType>& coeff,
    meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<diagStencil,Type,MeshType>> tSys
    (
        new linearSystem<diagStencil,Type,MeshType>(field)
    );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    Sys.A() = Zero;

    forAll(field, l)
    forAll(field[l], d)
    {
        meshDirection<diagStencil,MeshType>& A = Sys.A()[l][d];
        const meshDirection<scalar,MeshType>& c = coeff[l][d];
        const meshDirection<scalar,MeshType>& v = cv[l][d];

        forAllCells(A, i, j, k)
        {
            A(i,j,k) = v(i,j,k)*c(i,j,k);
        }
    }

    Sys.b() = Zero;

    return tSys;
}

}

}

}

}

#endif
