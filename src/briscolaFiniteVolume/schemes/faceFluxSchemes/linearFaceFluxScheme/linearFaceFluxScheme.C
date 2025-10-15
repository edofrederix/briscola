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

    const colocatedVectorFaceField& fan =
        this->fvMsh().metrics<colocated>().soa().faceAreaNormals();

    const colocatedScalarFaceField& fwc =
        this->fvMsh().metrics<colocated>().soa().faceWeightsCenter();

    const colocatedScalarFaceField& fwn =
        this->fvMsh().metrics<colocated>().soa().faceWeightsNeighbor();

    forAllFaces(Flux, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Flux(ijk)[fd*2] =
            (
                (
                    field(ijk)*fwc[fd](ijk)
                  + field(nei)*fwn[fd](ijk)
                )
              & fan[fd](ijk)
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

    const staggeredScalarFaceField& fa =
        this->fvMsh().metrics<staggered>().soa().faceAreas();

    const staggeredScalarFaceField& fwc =
        this->fvMsh().metrics<staggered>().soa().faceWeightsCenter();

    const staggeredScalarFaceField& fwn =
        this->fvMsh().metrics<staggered>().soa().faceWeightsNeighbor();

    forAllFaces(Flux, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        Flux(d,ijk)[fd*2] =
          - fa[fd](d,ijk)
          * (
                fwc[d](fd,ijk)*field(fd,ijk)
              + fwn[d](fd,ijk)*field(fd,lowerNeighbor(ijk,d))
            );

        Flux(d,nei)[fd*2+1] = -Flux(d,ijk)[fd*2];
    }

    return tFlux;
}

}

}

}
