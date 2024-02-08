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
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>
        (
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    A = Zero;
    b = Zero;

    forAllCells(A, l, d, i, j, k)
    {
        labelVector ijk(i,j,k);

        for (int fd = 0; fd < 3; fd++)
        {
            labelVector nei(upperNei(ijk,fd));

            const scalar lowerPhi =   factor*phi(l,d,ijk)[fd];
            const scalar upperPhi = - factor*phi(l,d,nei)[fd];

            A(l,d,ijk)[fd*2+1] = 0.5*lowerPhi;
            A(l,d,ijk)[fd*2+2] = 0.5*upperPhi;

            A(l,d,ijk)[0] += 0.5*(lowerPhi + upperPhi);
        }
    }

    const_cast<meshField<lowerFaceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::exDiv
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

    const meshField<scalar,MeshType>& cv =
        phi.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllCells(Div, d, i, j, k)
    {
        labelVector ijk(i,j,k);

        for (int fd = 0; fd < 3; fd++)
        {
            labelVector low(lowerNei(ijk,fd));
            labelVector upp(upperNei(ijk,fd));

            const scalar lowerPhi =   phi(d,ijk)[fd];
            const scalar upperPhi = - phi(d,upp)[fd];

            Div(d,ijk) +=
                (
                    0.5*(lowerPhi+upperPhi)*field(d,ijk)
                  + 0.5*lowerPhi*field(d,low)
                  + 0.5*upperPhi*field(d,upp)
                )
              / cv(d,ijk);
        }
    }

    return tDiv;
}

}

}

}
