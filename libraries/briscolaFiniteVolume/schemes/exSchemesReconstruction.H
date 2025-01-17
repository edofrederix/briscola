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
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Recon(i,j,k) = Zero;
        #endif

        for (int d = 0; d < 3; d++)
            Recon(i,j,k) +=
                0.5
              * (
                  - fn(i,j,k)[d*2  ]*field(d,i,j,k)
                  + fn(i,j,k)[d*2+1]*field(d,upperNeighbor(i,j,k,d))
                );
    }

    return tRecon;
}

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<typename outerProduct<vector,Type>::type,colocated>> tRecon
    (
        reconstruct(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tRecon;
}

// Reconstruct a centered field from a face flux

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
reconstruct
(
    const meshField<FaceSpace<Type>,MeshType>& flux
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

        for (label f = 0; f < 6; f++)
            weight += fn(d,i,j,k)[f]*fan(d,i,j,k)[f];

        ReconType sum = Zero;

        for (label f = 0; f < 6; f++)
            sum += fn(d,i,j,k)[f]*flux(d,i,j,k)[f];

        Recon(i,j,k) = inv(weight) & sum;
    }

    return tRecon;
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
reconstruct
(
    const tmp<meshField<FaceSpace<Type>,MeshType>>& tFlux
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,MeshType>> tRecon
    (
        reconstruct(tFlux())
    );

    if (tFlux.isTmp())
        tFlux.clear();

    return tRecon;
}

// Reconstruct a staggered field from a colocated face flux

template<class Type>
inline tmp<meshField<Type,staggered>> stagReconstruct
(
    const meshField<FaceSpace<Type>,colocated>& flux
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

    const meshField<faceScalar,colocated>& fa =
        flux.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(Recon, d, i, j, k)
        Recon(d,i,j,k) =
          - flux(i,j,k)[d*2]/fa(i,j,k)[d*2];

    return tRecon;
}

template<class Type>
inline tmp<meshField<Type,staggered>> stagReconstruct
(
    const tmp<meshField<FaceSpace<Type>,colocated>>& tFlux
)
{
    tmp<meshField<Type,staggered>> tRecon
    (
        stagReconstruct(tFlux())
    );

    if (tFlux.isTmp())
        tFlux.clear();

    return tRecon;
}

}

}

}

}

#endif
