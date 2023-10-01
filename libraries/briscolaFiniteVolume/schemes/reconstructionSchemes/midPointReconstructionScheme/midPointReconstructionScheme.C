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

    Recon = Zero;

    const meshField<faceVector,colocated>& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    forAllCells(Recon, i, j, k)
        Recon(i,j,k) =
            0.5
          * (
              - field(0,i,  j,  k  ) * fn(i,j,k).left()
              + field(0,i+1,j,  k  ) * fn(i,j,k).right()
              - field(1,i,  j,  k  ) * fn(i,j,k).bottom()
              + field(1,i,  j+1,k  ) * fn(i,j,k).top()
              - field(2,i,  j,  k  ) * fn(i,j,k).aft()
              + field(2,i,  j,  k+1) * fn(i,j,k).fore()
            );

    return tRecon;
}

}

}

}
