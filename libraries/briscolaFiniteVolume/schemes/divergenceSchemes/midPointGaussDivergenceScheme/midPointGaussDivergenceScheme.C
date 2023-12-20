#include "midPointGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointGaussDivergenceScheme<Type,MeshType>::midPointGaussDivergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    divergenceScheme<stencil,Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointGaussDivergenceScheme<Type,MeshType>::midPointGaussDivergenceScheme
(
    const fvMesh& fvMsh
)
:
    divergenceScheme<stencil,Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::imDiv
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

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar a = phi(l,d,ijk)[fd]*0.5;

        A(l,d,ijk).center() += a;
        A(l,d,nei).center() -= a;

        A(l,d,ijk)[1+fd*2  ] =   a;
        A(l,d,nei)[1+fd*2+1] = - a;
    }

    Sys.b() = Zero;

    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::exDiv
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

    const meshField<scalar,MeshType>& cv =
        phi.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllFaces(Div, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Type value = 0.5*(field(d,ijk) + field(d,nei));

        Div(d,ijk) += phi(d,ijk)[fd]*value/cv(d,ijk);
        Div(d,nei) -= phi(d,ijk)[fd]*value/cv(d,nei);
    }

    return tDiv;
}

}

}

}
