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
    const meshField<Type,staggered>& f
)
{
    typedef typename reconstructionScheme<Type>::ReconType ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+f.name()+")",
            f.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

    Recon = Zero;

    const meshField<faceVector,colocated>& fn =
        f.fvMsh().template metrics<colocated>().faceNormals();

    forAllLevels(Recon, l, d, i, j, k)
        Recon(l,d,i,j,k) =
            0.5
          * (
              - f(l,0,i,  j,  k  ) * fn(l,d,i,j,k).left()
              + f(l,0,i+1,j,  k  ) * fn(l,d,i,j,k).right()
              - f(l,1,i,  j,  k  ) * fn(l,d,i,j,k).bottom()
              + f(l,1,i,  j+1,k  ) * fn(l,d,i,j,k).top()
              - f(l,2,i,  j,  k  ) * fn(l,d,i,j,k).aft()
              + f(l,2,i,  j,  k+1) * fn(l,d,i,j,k).fore()
            );

    return tRecon;
}

}

}

}
