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
    const meshField<lowerFaceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    phi.restrict();

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

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllCells(A, l, d, i, j, k)
    {
        labelVector ijk(i,j,k);

        for (int fd = 0; fd < 3; fd++)
        {
            labelVector nei(upperNei(ijk,fd));

            const scalar lowerPhi =   factor*phi(l,d,ijk)[fd];
            const scalar upperPhi = - factor*phi(l,d,nei)[fd];

            A(l,d,ijk)[fd*2+1] = lowerPhi*fwn(l,d,ijk)[fd*2  ];
            A(l,d,ijk)[fd*2+2] = upperPhi*fwn(l,d,ijk)[fd*2+1];

            A(l,d,ijk)[0] +=
                lowerPhi*fwc(l,d,ijk)[fd*2  ]
              + upperPhi*fwc(l,d,ijk)[fd*2+1];
        }
    }

    phi.makeShallow();

    return tSys;
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
                    lowerPhi*fwc(d,ijk)[fd*2  ]*field(d,ijk)
                  + upperPhi*fwc(d,ijk)[fd*2+1]*field(d,ijk)
                  + lowerPhi*fwn(d,ijk)[fd*2  ]*field(d,low)
                  + upperPhi*fwn(d,ijk)[fd*2+1]*field(d,upp)
                )
              / cv(d,ijk);
        }
    }

    return tDiv;
}

}

}

}
