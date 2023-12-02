#include "linearFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

linearFaceFluxScheme::linearFaceFluxScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dict,fvMsh)
{}

linearFaceFluxScheme::linearFaceFluxScheme
(
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dictionary(),fvMsh)
{}

tmp<colocatedFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    tmp<colocatedFaceScalarField> tFlux
    (
        new colocatedFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    const meshField<faceScalar,colocated>& fwc =
        field.fvMsh().template metrics<colocated>().faceWeightsCenter();

    const meshField<faceScalar,colocated>& fwn =
        field.fvMsh().template metrics<colocated>().faceWeightsNeighbor();

    forAllCells(Flux, i, j, k)
        for (int f = 0; f < 6; f++)
            Flux(i,j,k)[f] =
                (
                    (
                        fwc(i,j,k)[f]*field(i,j,k)
                      + fwn(i,j,k)[f]*field(labelVector(i,j,k)+faceOffsets[f])
                    )
                  & fan(i,j,k)[f]
                );

    return tFlux;
}

tmp<staggeredFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredFaceScalarField> tFlux
    (
        new staggeredFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    staggeredFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    const meshField<faceScalar,staggered>& fwc =
        field.fvMsh().template metrics<staggered>().faceWeightsCenter();

    const meshField<faceScalar,staggered>& fwn =
        field.fvMsh().template metrics<staggered>().faceWeightsNeighbor();

    forAllDirections(Flux, d, i, j, k)
    {
        for (int f = 0; f < 6; f++)
        {
            const labelVector ijk(i,j,k);

            const label e = f/2;
            const label g = (e == d ? f : d*2);

            const labelVector c(e == d || f%2 == 0 ? ijk : ijk+faceOffsets[f]);
            const labelVector n(c + faceOffsets[g]);

            const label sign = (f%2)*2-1;

            Flux(d,ijk)[f] =
                fa(d,ijk)[f]
              * sign
              * (
                    fwc(e,c)[g]*field(e,c)
                  + fwn(e,c)[g]*field(e,n)
                );
        }
    }

    return tFlux;
}

}

}

}
