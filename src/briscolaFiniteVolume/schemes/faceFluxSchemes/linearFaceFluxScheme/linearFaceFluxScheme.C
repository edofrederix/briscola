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
    faceField
    <
        typename innerProduct<Type,vector>::type,
        colocated
    >
> linearFaceFluxScheme::coloFaceFlux
(
    const meshField<Type,colocated>& field
)
{
    typedef
        faceField
        <
            typename innerProduct<Type,vector>::type,
            colocated
        > returnType;

    tmp<returnType> tPhi =
        returnType::New
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        );

    returnType& phi = tPhi.ref();

    const colocatedVectorFaceField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    const colocatedScalarFaceField& fwc =
        this->fvMsh().metrics<colocated>().faceWeightsCenter();

    const colocatedScalarFaceField& fwn =
        this->fvMsh().metrics<colocated>().faceWeightsNeighbor();

    forAllFaces(phi, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        phi[fd](ijk) =
            (
                (
                    field(ijk)*fwc[fd](ijk)
                  + field(nei)*fwn[fd](ijk)
                )
              & fan[fd](ijk)
            );
    }

    return tPhi;
}

linearFaceFluxScheme::linearFaceFluxScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    faceFluxScheme(fvMsh, is)
{}

tmp<colocatedScalarFaceField> linearFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<colocatedVectorFaceField> linearFaceFluxScheme::faceFlux
(
    const colocatedTensorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<staggeredScalarFaceField> linearFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredScalarFaceField> tPhi =
        staggeredScalarFaceField::New
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        );

    staggeredScalarFaceField& phi = tPhi.ref();

    const staggeredScalarFaceField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    const staggeredScalarFaceField& fwc =
        this->fvMsh().metrics<staggered>().faceWeightsCenter();

    const staggeredScalarFaceField& fwn =
        this->fvMsh().metrics<staggered>().faceWeightsNeighbor();

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        phi[fd](d,ijk) =
          - fa[fd](d,ijk)
          * (
                fwc[d](fd,ijk)*field(fd,ijk)
              + fwn[d](fd,ijk)*field(fd,nei)
            );
    }

    return tPhi;
}

}

}

}
