#include "midPointGaussDivergenceScheme.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointGaussDivergenceScheme<Type,MeshType>::midPointGaussDivergenceScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    divergenceScheme<stencil,Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::imDiv
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

    meshField<stencil,MeshType>& A = sys.A();

    A = Zero;

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const scalar value = 0.5*phi[fd](l,d,ijk);

        A(l,d,ijk)[fd*2+1] =  value;
        A(l,d,nei)[fd*2+2] = -value;

        A(l,d,ijk)[0] += value;
        A(l,d,nei)[0] -= value;
    }

    sys.b() = Zero;

    if (shallow)
        collapse(phi);

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::exDiv
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

    const meshField<scalar,MeshType>& icv =
        phi.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    Div = Zero;

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value =
            phi[fd](d,ijk)*0.5*(field(d,ijk) + field(d,nei));

        Div(d,ijk) += value;
        Div(d,nei) -= value;
    }

    Div *= icv;

    return tDiv;
}

}

}

}
