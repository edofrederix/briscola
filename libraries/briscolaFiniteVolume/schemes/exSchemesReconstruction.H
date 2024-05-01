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

    Recon = Zero;

    const meshField<faceVector,colocated>& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    forAllCells(Recon, i, j, k)
        for (int d = 0; d < 3; d++)
            Recon(i,j,k) +=
                0.5
              * (
                  - fn(i,j,k)[d*2  ]*field(d,i,j,k)
                  + fn(i,j,k)[d*2+1]*field(d,upperNei(i,j,k,d))
                );

    return tRecon;
}

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const tmp<meshField<Type,staggered>>& tfield
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        reconstruct(tfield())
    );

    tfield.clear();

    return tRecon;
}

// Reconstruct a colocated field from a face field

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<LowerFaceSpace<Type>,colocated>& field
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

    const meshField<faceVector,colocated>& fan =
        field.fvMsh().template metrics<colocated>().faceAreaNormals();

    forAllCells(Recon, i, j, k)
    {
        tensor weight = Zero;

        for (int f = 0; f < 6; f++)
            weight += fn(i,j,k)[f]*fan(i,j,k)[f];

        ReconType sum = Zero;

        for (int d = 0; d < 3; d++)
            sum +=
                fn(i,j,k)[d*2  ]*field(i,j,k)[d]
              - fn(i,j,k)[d*2+1]*field(upperNei(i,j,k,d))[d];

        Recon(i,j,k) = inv(weight) & sum;
    }

    return tRecon;
}

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const tmp<meshField<LowerFaceSpace<Type>,colocated>>& tfield
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        reconstruct(tfield())
    );

    tfield.clear();

    return tRecon;
}

}

}

}

}

#endif
