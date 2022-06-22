#include "midPointFaceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointFaceGradientScheme<Type,MeshType>::midPointFaceGradientScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointFaceGradientScheme<Type,MeshType>::midPointFaceGradientScheme
(
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<Hex<Type>,MeshType>>
midPointFaceGradientScheme<Type,MeshType>::faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Hex<Type>,MeshType>> tGrad
    (
        new meshField<Hex<Type>,MeshType>
        (
            "faceGrad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Hex<Type>,MeshType>& Grad = tGrad.ref();

    forAll(field, l)
    forAll(field[l], d)
    {
        const meshDirection<hexScalar,MeshType>& fd =
            field.fvMsh().template
            metrics<MeshType>().faceDeltas()[l][d];

        meshDirection<Hex<Type>,MeshType>& G = Grad[l][d];
        const meshDirection<Type,MeshType>& f = field[l][d];

        G.initGhosts();

        forAllCells(f, i, j, k)
        {
            G(i,j,k) =
                Hex<Type>
                (
                    f(i,j,k)-f(i-1,j,  k  ),
                  - f(i,j,k)+f(i+1,j,  k  ),
                    f(i,j,k)-f(i,  j-1,k  ),
                  - f(i,j,k)+f(i,  j+1,k  ),
                    f(i,j,k)-f(i,  j,  k-1),
                  - f(i,j,k)+f(i,  j,  k+1)
                )
              * fd(i,j,k);
        }
    }

    return tGrad;
}

}

}

}
