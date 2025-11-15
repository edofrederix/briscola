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
    const bool shallowPhi = phi.shallow();
    const bool shallowField = field.shallow();

    if (shallowPhi)
        restrict(phi);

    if (shallowField)
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

    faceField<scalar,MeshType> psi(limiter_->psi(phi,field,true));

    A = Zero;

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const scalar w1 = psi[fd](l,d,ijk);
        const scalar w2 = (1.0 - w1);

        const scalar flux = phi[fd](l,d,ijk);

        const scalar value1 =
            flux*(w1*fwn[fd](l,d,ijk) + w2*(flux < 0));

        const scalar value2 =
            flux*(w1*fwc[fd](l,d,ijk) + w2*(flux >= 0));

        A(l,d,ijk)[fd*2+1] =  value1;
        A(l,d,nei)[fd*2+2] = -value2;

        A(l,d,ijk)[0] += value2;
        A(l,d,nei)[0] -= value1;
    }

    sys.b() = Zero;

    if (shallowPhi)
        collapse(phi);

    if (shallowField)
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

    const faceField<scalar,MeshType> psi(limiter_->psi(phi,field,false));

    Div = Zero;

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const scalar w1 = psi[fd](d,ijk);
        const scalar w2 = (1.0 - w1);

        const scalar flux = phi[fd](d,ijk);

        const Type value =
            flux
          * (
                field(d,ijk)*(w1*fwc[fd](d,ijk) + w2*(flux >= 0))
              + field(d,nei)*(w1*fwn[fd](d,ijk) + w2*(flux < 0))
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
