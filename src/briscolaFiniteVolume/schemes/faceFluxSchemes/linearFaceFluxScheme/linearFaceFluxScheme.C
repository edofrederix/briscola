#include "linearFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
tmp
<
    meshField
    <
        FaceSpace<typename innerProduct<Type,vector>::type>,
        colocated
    >
> linearFaceFluxScheme::coloFaceFlux
(
    const meshField<Type,colocated>& field
)
{
    typedef
        meshField
        <
            FaceSpace<typename innerProduct<Type,vector>::type>,
            colocated
        > returnType;

    tmp<returnType> tFlux =
        returnType::New
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        );

    returnType& Flux = tFlux.ref();

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    const colocatedFaceScalarField& fwc =
        field.fvMsh().template metrics<colocated>().faceWeightsCenter();

    const colocatedFaceScalarField& fwn =
        field.fvMsh().template metrics<colocated>().faceWeightsNeighbor();

    forAllFaces(Flux, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Flux(ijk)[fd*2] =
            (
                (
                    field(ijk)*fwc(ijk)[fd*2]
                  + field(nei)*fwn(ijk)[fd*2]
                )
              & fan(ijk)[fd*2]
            );

        Flux(nei)[fd*2+1] = -Flux(ijk)[fd*2];
    }

    return tFlux;
}

linearFaceFluxScheme::linearFaceFluxScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    faceFluxScheme(fvMsh, is)
{}

tmp<colocatedFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<colocatedFaceVectorField> linearFaceFluxScheme::faceFlux
(
    const colocatedTensorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<staggeredFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredFaceScalarField> tFlux =
        staggeredFaceScalarField::New
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        );

    staggeredFaceScalarField& Flux = tFlux.ref();

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    const staggeredFaceScalarField& fwc =
        field.fvMsh().template metrics<staggered>().faceWeightsCenter();

    const staggeredFaceScalarField& fwn =
        field.fvMsh().template metrics<staggered>().faceWeightsNeighbor();

    forAllFaces(Flux, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        Flux(d,ijk)[fd*2] =
          - fa(d,ijk)[fd*2]
          * (
                fwc(fd,ijk)[d*2]*field(fd,ijk)
              + fwn(fd,ijk)[d*2]*field(fd,lowerNeighbor(ijk,d))
            );

        Flux(d,nei)[fd*2+1] = -Flux(d,ijk)[fd*2];
    }

    return tFlux;
}

}

}

}
