#include "midPointFaceFluxScheme.H"

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
> midPointFaceFluxScheme::coloFaceFlux
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

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        phi[fd](l,d,ijk) =
            0.5*((field(l,d,ijk) + field(l,d,nei)) & fan[fd](l,d,ijk));
    }

    return tPhi;
}

midPointFaceFluxScheme::midPointFaceFluxScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    faceFluxScheme(fvMsh, is)
{}

tmp<colocatedScalarFaceField> midPointFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<colocatedVectorFaceField> midPointFaceFluxScheme::faceFlux
(
    const colocatedTensorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<staggeredScalarFaceField> midPointFaceFluxScheme::faceFlux
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

    forAllFaces(phi, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        phi[fd](l,d,ijk) =
          - 0.5*fa[fd](l,d,ijk)*(field(l,fd,ijk) + field(l,fd,nei));
    }

    return tPhi;
}

}

}

}
