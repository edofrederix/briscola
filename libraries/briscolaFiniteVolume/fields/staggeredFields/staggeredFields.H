#ifndef staggeredFields_H
#define staggeredFields_H

#include "meshFields.H"
#include "staggered.H"
#include "staggeredFieldsFwd.H"
#include "staggeredLevels.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Special function to perform a dot product between a vector dimension (first
// argument) and a staggered dimension (second argument). The vector must be
// given at staggered locations. When the staggered dimension index is not equal
// to the vector dimension index, interpolation is required.

inline tmp<staggeredScalarField>
stagDotProduct
(
    const staggeredVectorField& f,
    const staggeredVectorField& g
)
{
    tmp<staggeredScalarField> tProd
    (
        new staggeredScalarField
        (
            "("+f.name()+"&"+g.name()+")",
            f.fvMsh()
        )
    );

    staggeredScalarField& Prod = tProd.ref();

    Prod = Zero;

    forAllCells(Prod, di, i, j, k)
        for (int dj = 0; dj < 3; dj++)
            Prod(di,i,j,k) +=
                f(di,i,j,k)[dj]
              * (
                    di == dj
                  ? g(dj,i,j,k)[di]
                  : 0.25
                  * (
                        g(dj,i,j,k)[di]
                      + g(dj,lowerNei(i,j,k,di))[di]
                      + g(dj,upperNei(i,j,k,dj))[di]
                      + g(dj,lowerNei(upperNei(i,j,k,dj),di))[di]
                    )
                );

    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const tmp<staggeredVectorField>& tf,
    const staggeredVectorField& g
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(tf(),g));
    tf.clear();
    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const staggeredVectorField& f,
    const tmp<staggeredVectorField>& tg
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(f,tg()));
    tg.clear();
    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const tmp<staggeredVectorField>& tf,
    const tmp<staggeredVectorField>& tg
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(tf(),tg()));
    tf.clear();
    tg.clear();
    return tProd;
}

}

}

}

#endif
