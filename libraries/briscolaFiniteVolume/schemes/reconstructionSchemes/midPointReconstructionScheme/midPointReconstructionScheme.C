#include "midPointReconstructionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
midPointReconstructionScheme<Type>::midPointReconstructionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    reconstructionScheme<Type>(dict,fvMsh)
{}

template<class Type>
midPointReconstructionScheme<Type>::midPointReconstructionScheme
(
    const fvMesh& fvMsh
)
:
    reconstructionScheme<Type>(dictionary(),fvMsh)
{}

template<class Type>
tmp<meshField<typename reconstructionScheme<Type>::ReconType, colocated>>
midPointReconstructionScheme<Type>::reconstruct
(
    const meshField<Type,staggered>& field
)
{
    typedef typename reconstructionScheme<Type>::ReconType ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

    forAll(field, l)
    {
        const meshDirection<faceVector,colocated>& fn =
            field.fvMsh().template
            metrics<colocated>().faceNormals()[l][0];

        meshDirection<ReconType,colocated>& R = Recon[l][0];

        const meshDirection<Type,staggered>& f0 = field[l][0];
        const meshDirection<Type,staggered>& f1 = field[l][1];
        const meshDirection<Type,staggered>& f2 = field[l][2];

        R.initGhosts();

        forAllCells(R, i, j, k)
        {
            R(i,j,k) =
                0.5
              * (
                    f0(i,  j,  k  )*fn(i,j,k).left()
                  + f0(i+1,j,  k  )*fn(i,j,k).right()
                  + f1(i,  j,  k  )*fn(i,j,k).bottom()
                  + f1(i,  j+1,k  )*fn(i,j,k).top()
                  + f2(i,  j,  k  )*fn(i,j,k).aft()
                  + f2(i,  j,  k+1)*fn(i,j,k).fore()
                );
        }
    }

    return tRecon;
}

}

}

}
