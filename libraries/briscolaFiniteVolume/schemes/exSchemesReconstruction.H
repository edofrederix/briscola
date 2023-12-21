#ifndef exSchemesReconstruct_H
#define exSchemesReconstruct_H

#include "meshField.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Reconstruct a colocated field from a staggered one

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<Type,staggered>& field
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

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

}

#endif
