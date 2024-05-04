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
    const tmp<meshField<Type,staggered>>& tField
)
{
    tmp<meshField<typename outerProduct<vector,Type>::type,colocated>> tRecon
    (
        reconstruct(tField())
    );

    tField.clear();

    return tRecon;
}

// Reconstruct a centered field from a face flux

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
reconstruct
(
    const meshField<LowerFaceSpace<Type>,MeshType>& flux
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,MeshType>> tRecon
    (
        new meshField<ReconType,MeshType>
        (
            "reconstruct("+flux.name()+")",
            flux.fvMsh()
        )
    );

    meshField<ReconType,MeshType>& Recon = tRecon.ref();

    const meshField<faceVector,MeshType>& fn =
        flux.fvMsh().template metrics<MeshType>().faceNormals();

    const meshField<faceVector,MeshType>& fan =
        flux.fvMsh().template metrics<MeshType>().faceAreaNormals();

    forAllCells(Recon, d, i, j, k)
    {
        tensor weight = Zero;

        for (int f = 0; f < 6; f++)
            weight += fn(d,i,j,k)[f]*fan(d,i,j,k)[f];

        ReconType sum = Zero;

        for (int fd = 0; fd < 3; fd++)
            sum +=
                fn(d,i,j,k)[fd*2  ]*flux(d,i,j,k)[fd]
              - fn(d,i,j,k)[fd*2+1]*flux(d,upperNei(i,j,k,fd))[fd];

        Recon(i,j,k) = inv(weight) & sum;
    }

    return tRecon;
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
reconstruct
(
    const tmp<meshField<LowerFaceSpace<Type>,MeshType>>& tFlux
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,MeshType>> tRecon
    (
        reconstruct(tFlux())
    );

    tFlux.clear();

    return tRecon;
}

// Reconstruct a staggered field from a colocated face flux

template<class Type>
inline tmp<meshField<Type,staggered>> stagReconstruct
(
    const meshField<LowerFaceSpace<Type>,colocated>& flux
)
{
    tmp<meshField<Type,staggered>> tRecon
    (
        new meshField<Type,staggered>
        (
            "reconstruct("+flux.name()+")",
            flux.fvMsh()
        )
    );

    meshField<Type,staggered>& Recon = tRecon.ref();

    Recon = Zero;

    const meshField<faceScalar,colocated>& fa =
        flux.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(Recon, d, i, j, k)
        Recon(d,i,j,k) =
          - flux(i,j,k)[d]/fa(i,j,k)[d*2];

    return tRecon;
}

template<class Type>
inline tmp<meshField<Type,staggered>> stagReconstruct
(
    const tmp<meshField<LowerFaceSpace<Type>,colocated>>& tFlux
)
{
    tmp<meshField<Type,staggered>> tRecon
    (
        stagReconstruct(tFlux())
    );

    tFlux.clear();

    return tRecon;
}

}

}

}

}

#endif
