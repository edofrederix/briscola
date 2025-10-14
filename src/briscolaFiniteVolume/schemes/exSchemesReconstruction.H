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

    tmp<meshField<ReconType,colocated>> tRecon =
        meshField<ReconType,colocated>::New
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
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
              - fn(i,j,k)[d*2  ]*field(d,i,j,k)
              + fn(i,j,k)[d*2+1]*field(d,upperNeighbor(i,j,k,d));

        Recon(i,j,k) *= 0.5;
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

    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon = reconstruct(tField());

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
    typedef meshField<tensor,MeshType> WeightFieldType;

    const fvMesh& fvMsh = flux.fvMsh();

    tmp<meshField<ReconType,MeshType>> tRecon =
        meshField<ReconType,MeshType>::New
        (
            "reconstruct("+flux.name()+")",
            flux.fvMsh()
        );

    meshField<ReconType,MeshType>& Recon = tRecon.ref();

    const meshField<faceVector,MeshType>& fn =
        fvMsh.template metrics<MeshType>().faceNormals();

    // Create reconstruction weights and store

    if (!fvMsh.db().foundObject<WeightFieldType>("reconstructionWeights"))
    {
        // Weight field

        tmp<WeightFieldType> tWeights =
            WeightFieldType::New
            (
                "reconstructionWeights",
                fvMsh
            );

        WeightFieldType& weights = tWeights.ref();

        const meshField<faceVector,MeshType>& fan =
            flux.fvMsh().template metrics<MeshType>().faceAreaNormals();

        // Compute weight matrix and invert

        forAllCells(weights, d, i, j, k)
        {
            #ifdef NO_BLOCK_ZERO_INIT
            weights(d,i,j,k) = Zero;
            #endif

            for (label f = 0; f < 6; f++)
                weights(d,i,j,k) +=
                    fn(d,i,j,k)[f]*fan(d,i,j,k)[f];

            weights(d,i,j,k) = inv(weights(d,i,j,k));
        }

        // Store

        fvMsh.db().objectRegistry::store(tWeights.ptr());
    }

    const WeightFieldType& weights =
        fvMsh.db().lookupObject<WeightFieldType>("reconstructionWeights");

    forAllCells(Recon, d, i, j, k)
    {
        ReconType sum = Zero;

        for (label f = 0; f < 6; f++)
            sum += fn(d,i,j,k)[f]*flux(d,i,j,k)[f];

        Recon(d,i,j,k) = weights(d,i,j,k) & sum;
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

    tmp<meshField<ReconType,MeshType>> tRecon = reconstruct(tFlux());

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
    tmp<meshField<Type,staggered>> tRecon =
        meshField<Type,staggered>::New
        (
            "reconstruct("+flux.name()+")",
            flux.fvMsh()
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
    tmp<meshField<Type,staggered>> tRecon =
        stagReconstruct(tFlux());

    if (tFlux.isTmp())
        tFlux.clear();

    return tRecon;
}

}

}

}

}

#endif
