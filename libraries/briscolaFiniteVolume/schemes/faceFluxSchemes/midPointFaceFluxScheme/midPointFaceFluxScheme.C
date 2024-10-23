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
    meshField
    <
        LowerFaceSpace<typename innerProduct<Type,vector>::type>,
        colocated
    >
> midPointFaceFluxScheme::coloFaceFlux
(
    const meshField<Type,colocated>& field
)
{
    typedef
        meshField
        <
            LowerFaceSpace<typename innerProduct<Type,vector>::type>,
            colocated
        > returnType;

    tmp<returnType> tFlux
    (
        new returnType
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    returnType& Flux = tFlux.ref();

    Flux = Zero;

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    forAllFaces(Flux, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(lowerNei(i,j,k,fd));

        Flux(ijk)[fd] = (0.5*(field(ijk) + field(nei)) & fan(ijk)[fd*2]);
    }

    return tFlux;
}

midPointFaceFluxScheme::midPointFaceFluxScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    faceFluxScheme(fvMsh, is)
{}

tmp<colocatedLowerFaceScalarField> midPointFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<colocatedLowerFaceVectorField> midPointFaceFluxScheme::faceFlux
(
    const colocatedTensorField& field
)
{
    return this->coloFaceFlux(field);
}

tmp<staggeredLowerFaceScalarField> midPointFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredLowerFaceScalarField> tFlux
    (
        new staggeredLowerFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    staggeredLowerFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    forAllFaces(Flux, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(lowerNei(i,j,k,d));

        Flux(d,ijk)[fd] =
          - 0.5*fa(d,ijk)[fd*2] * (field(fd,ijk) + field(fd,nei));
    }

    return tFlux;
}

}

}

}
