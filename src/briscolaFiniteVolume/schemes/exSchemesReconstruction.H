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

    const colocatedVectorFaceField& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    Recon = Zero;

    forAllFaces(fn, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const ReconType value = -0.5*fn[fd](ijk)*field(fd,ijk);

        Recon(ijk) += value;
        Recon(nei) += value;
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
    const faceField<Type,MeshType>& flux
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

    meshField<ReconType,MeshType>& recon = tRecon.ref();

    const faceField<vector,MeshType>& fn =
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

        const faceField<vector,MeshType>& fan =
            flux.fvMsh().template metrics<MeshType>().faceAreaNormals();

        weights = Zero;

        forAllFaces(fn, fd, d, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(i,j,k,fd));

            const tensor value = fn[fd](d,ijk)*fan[fd](d,ijk);

            weights(d,ijk) += value;
            weights(d,nei) += value;
        }

        forAllCells(weights, d, i, j, k)
            weights(d,i,j,k) = inv(weights(d,i,j,k));

        // Store

        fvMsh.db().objectRegistry::store(tWeights.ptr());
    }

    const WeightFieldType& weights =
        fvMsh.db().lookupObject<WeightFieldType>("reconstructionWeights");

    recon = Zero;

    forAllFaces(fn, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const ReconType value = fn[fd](d,ijk)*flux[fd](d,ijk);

        recon(d,ijk) += value;
        recon(d,nei) += value;
    }

    forAllCells(recon, d, i, j, k)
        recon(d,i,j,k) = (weights(d,i,j,k) & recon(d,i,j,k));

    return tRecon;
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
reconstruct
(
    const tmp<faceField<Type,MeshType>>& tFlux
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
    const faceField<Type,colocated>& flux
)
{
    tmp<meshField<Type,staggered>> tRecon =
        meshField<Type,staggered>::New
        (
            "reconstruct("+flux.name()+")",
            flux.fvMsh()
        );

    meshField<Type,staggered>& recon = tRecon.ref();

    const faceField<scalar,colocated>& fa =
        flux.fvMsh().template metrics<colocated>().faceAreas();

    // Important in order not to set invalid staggered ghost cells
    recon = Zero;

    forAllCells(recon, d, i, j, k)
        recon(d,i,j,k) = -flux[d](i,j,k)/fa[d](i,j,k);

    return tRecon;
}

template<class Type>
inline tmp<meshField<Type,staggered>> stagReconstruct
(
    const tmp<faceField<Type,colocated>>& tFlux
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
