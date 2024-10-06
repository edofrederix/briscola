#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
linearGaussLaplacianScheme<SType,Type,MeshType>::linearGaussLaplacianScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    laplacianScheme<SType,Type,MeshType>(fvMsh, is)
{}

template<class SType, class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussLaplacianScheme<SType,Type,MeshType>::exLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    tmp<meshField<Type,MeshType>> tLap
    (
        new meshField<Type,MeshType>
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,MeshType>& Lap = tLap.ref();

    Lap = Zero;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllCells(Lap, d, i, j, k)
    {
        labelVector ijk(i,j,k);

        for (int fd = 0; fd < 3; fd++)
        {
            labelVector low(lowerNei(ijk,fd));
            labelVector upp(upperNei(ijk,fd));

            scalar lower = factor*fa(d,ijk)[fd*2  ]*delta(d,ijk)[fd*2  ];
            scalar upper = factor*fa(d,ijk)[fd*2+1]*delta(d,ijk)[fd*2+1];

            if (lambdaPtr)
            {
                lower *= lambdaPtr->operator()(d,ijk)[fd];
                upper *= lambdaPtr->operator()(d,upp)[fd];
            }

            Lap(d,ijk) +=
                (
                    lower*(field(d,low) - field(d,ijk))
                  + upper*(field(d,upp) - field(d,ijk))
                )
              / cv(d,ijk);
        }
    }

    return tLap;
}

}

}

}
