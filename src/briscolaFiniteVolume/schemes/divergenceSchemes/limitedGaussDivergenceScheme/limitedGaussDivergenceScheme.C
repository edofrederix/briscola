#include "limitedGaussDivergenceScheme.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
limitedGaussDivergenceScheme<Type,MeshType>::limitedGaussDivergenceScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    divergenceScheme<stencil,Type,MeshType>(fvMsh, is),
    limiter_(limiterScheme<Type,MeshType>::New(fvMsh,is))
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
limitedGaussDivergenceScheme<Type,MeshType>::imDiv
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    restrict(phi);
    restrict(field);

    tmp<linearSystem<stencil,Type,MeshType>> tSys =
        linearSystem<stencil,Type,MeshType>::New
        (
            word("div(" + phi.name() + "," + field.name() + ")"),
            const_cast<meshField<Type,MeshType>&>(field)
        );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();

    const faceField<scalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const faceField<scalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const faceField<scalar,MeshType> psi(limiter_->psi(phi,field));

    restrict(psi);

    A = Zero;

    forAllCells(A, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        for (label f = 0; f < 6; f++)
        {
            const labelVector upp(upperFaceNeighbor(ijk,f));
            const label fd = f/2;

            A(l,d,ijk)[f+1] =
                phi[fd](l,d,upp)
              * (
                    psi[fd](l,d,upp)*fwn[fd](l,d,upp)
                  + (1.0 - psi[fd](l,d,upp))*(phi[fd](l,d,upp) < 0)
                );

            A(l,d,ijk)[0] +=
                phi[fd](l,d,upp)
              * (
                    psi[fd](l,d,upp)*fwc[fd](l,d,upp)
                  + (1.0 - psi[fd](l,d,upp))*(phi[fd](l,d,upp) >= 0)
                );
        }
    }

    sys.b() = Zero;

    collapse(phi);
    collapse(field);

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
limitedGaussDivergenceScheme<Type,MeshType>::exDiv
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

    const faceField<scalar,MeshType> psi(limiter_->psi(phi,field));

    Div = Zero;

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value =
            phi[fd](d,ijk)
          * (
                field(d,ijk)
              * (
                    psi[fd](d,ijk)*fwc[fd](d,ijk)
                  + (1.0 - psi[fd](d,ijk))*(phi[fd](d,ijk) >= 0)
                )
              + field(d,nei)
              * (
                    psi[fd](d,ijk)*fwn[fd](d,ijk)
                  + (1.0 - psi[fd](d,ijk))*(phi[fd](d,ijk) < 0)
                )
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
