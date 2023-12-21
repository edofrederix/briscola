#include "linearGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussDivergenceScheme<Type,MeshType>::linearGaussDivergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    divergenceScheme<stencil,Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearGaussDivergenceScheme<Type,MeshType>::linearGaussDivergenceScheme
(
    const fvMesh& fvMsh
)
:
    divergenceScheme<stencil,Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
template<template<class> class OpType>
void linearGaussDivergenceScheme<Type,MeshType>::imDiv
(
    linearSystem<stencil,Type,MeshType>& sys,
    const meshField<lowerFaceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).restrict();

    meshField<stencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    if (isEqOp<OpType>() || isEqMinusOp<OpType>())
    {
        A = Zero;
        b = Zero;
    }

    OpType<scalar> bop;

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar a = factor*phi(l,d,ijk)[fd]*fwc(l,d,ijk)[fd*2];
        scalar b = factor*phi(l,d,ijk)[fd]*fwn(l,d,ijk)[fd*2];

        bop(A(l,d,ijk)[1+fd*2  ],  b);
        bop(A(l,d,nei)[1+fd*2+1], -a);

        if (isEqOp<OpType>())
        {
            A(l,d,ijk)[0] += a;
            A(l,d,nei)[0] -= b;
        }
        else if (isEqMinusOp<OpType>())
        {
            A(l,d,ijk)[0] -= a;
            A(l,d,nei)[0] += b;
        }
        else
        {
            bop(A(l,d,ijk)[0],  a);
            bop(A(l,d,nei)[0], -b);
        }
    }

    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).makeShallow();
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::exDiv
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tDiv
    (
        new meshField<Type,MeshType>
        (
            "div("+phi.name()+","+field.name()+")",
            phi.fvMsh()
        )
    );

    meshField<Type,MeshType>& Div = tDiv.ref();

    Div = Zero;

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<scalar,MeshType>& cv =
        phi.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllFaces(Div, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Type value =
            fwc(d,ijk)[fd*2]*field(d,ijk)
          + fwn(d,ijk)[fd*2]*field(d,nei);

        Div(d,ijk) += phi(d,ijk)[fd]*value/cv(d,ijk);
        Div(d,nei) -= phi(d,ijk)[fd]*value/cv(d,nei);
    }

    return tDiv;
}

}

}

}
