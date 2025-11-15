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

    // We need to handle deep fields properly
    phi.make(field.deep());

    const colocatedVectorFaceField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    const colocatedScalarFaceField& fwc =
        this->fvMsh().metrics<colocated>().faceWeightsCenter();

    const colocatedScalarFaceField& fwn =
        this->fvMsh().metrics<colocated>().faceWeightsNeighbor();

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        phi[fd](l,d,ijk) =
            (
                (
                    field(l,d,ijk)*fwc[fd](l,d,ijk)
                  + field(l,d,nei)*fwn[fd](l,d,ijk)
                )
              & fan[fd](l,d,ijk)
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

    // We need to handle deep fields properly
    phi.make(field.deep());

    const staggeredScalarFaceField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    const staggeredScalarFaceField& fwc =
        this->fvMsh().metrics<staggered>().faceWeightsCenter();

    const staggeredScalarFaceField& fwn =
        this->fvMsh().metrics<staggered>().faceWeightsNeighbor();

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        phi[fd](l,d,ijk) =
          - fa[fd](l,d,ijk)
          * (
                fwc[d](l,fd,ijk)*field(l,fd,ijk)
              + fwn[d](l,fd,ijk)*field(l,fd,nei)
            );
    }

    return tPhi;
}

}

}

}
