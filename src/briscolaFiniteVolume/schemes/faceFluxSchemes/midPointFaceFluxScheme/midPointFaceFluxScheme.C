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

    const colocatedVectorFaceField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    forAllFaces(phi, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        phi[fd](ijk) = 0.5*((field(ijk) + field(nei)) & fan[fd](ijk));
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

    const staggeredScalarFaceField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        // Fluxes in x-direction should receive from the first field direction,
        // fluxes in y-direction from the second, etc

        phi[fd](d,ijk) =
          - 0.5*fa[fd](d,ijk)
          * (
                field(fd,ijk)
              + field(fd,lowerNeighbor(ijk,d))
            );
    }

    return tPhi;
}

}

}

}
