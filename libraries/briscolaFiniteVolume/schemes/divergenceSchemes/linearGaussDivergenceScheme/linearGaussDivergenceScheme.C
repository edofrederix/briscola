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
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::imDiv
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    meshField<stencil,MeshType>& A = Sys.A();
    A = Zero;

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar a = phi(l,d,ijk)[fd]*fwc(l,d,ijk)[fd*2];
        scalar b = phi(l,d,ijk)[fd]*fwn(l,d,ijk)[fd*2];

        A(l,d,ijk).center() += a;
        A(l,d,nei).center() -= b;

        A(l,d,ijk)[1+fd*2  ] =   b;
        A(l,d,nei)[1+fd*2+1] = - a;
    }

    Sys.b() = Zero;

    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::exDiv
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
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
