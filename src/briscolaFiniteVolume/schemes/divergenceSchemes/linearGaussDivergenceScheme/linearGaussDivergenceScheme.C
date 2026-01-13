#include "linearGaussDivergenceScheme.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussDivergenceScheme<Type,MeshType>::linearGaussDivergenceScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    divergenceScheme<stencil,Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::imDiv
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    const bool shallow = phi.shallow();

    if (shallow)
        restrict(phi);

    tmp<linearSystem<stencil,Type,MeshType>> tSys =
        linearSystem<stencil,Type,MeshType>::New
        (
            word("div(" + phi.name() + "," + field.name() + ")"),
            const_cast<meshField<Type,MeshType>&>(field)
        );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    sys.symmetric() = false;
    sys.diagonal() = false;

    meshField<stencil,MeshType>& A = sys.A();

    const faceField<scalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const faceField<scalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    A = Zero;

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const scalar value1 = phi[fd](l,d,ijk)*fwn[fd](l,d,ijk);
        const scalar value2 = phi[fd](l,d,ijk)*fwc[fd](l,d,ijk);

        A(l,d,ijk)[fd*2+1] =  value1;
        A(l,d,nei)[fd*2+2] = -value2;

        A(l,d,ijk)[0] += value2;
        A(l,d,nei)[0] -= value1;
    }

    sys.b() = Zero;

    if (shallow)
        collapse(phi);

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::exDiv
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tDiv =
        meshField<Type,MeshType>::New
        (
            "div("+phi.name()+","+field.name()+")",
            phi.fvMsh()
        );

    meshField<Type,MeshType>& Div = tDiv.ref();

    const faceField<scalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const faceField<scalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<scalar,MeshType>& icv =
        phi.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    Div = Zero;

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value =
            phi[fd](d,ijk)
          * (
                fwc[fd](d,ijk)*field(d,ijk)
              + fwn[fd](d,ijk)*field(d,nei)
            );

        Div(d,ijk) += value;
        Div(d,nei) -= value;
    }

    Div *= icv;

    return tDiv;
}

}

}

}
