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

// Special function to perform a dot product between a staggered dimension and
// vector dimension. When the staggered dimension index is not equal to the
// vector dimension index, interpolation is required.

inline tmp<staggeredScalarField>
stagDotProduct
(
    const staggeredVectorField& f1,
    const staggeredVectorField& f2
)
{
    tmp<staggeredScalarField> tProd
    (
        new staggeredScalarField
        (
            "("+f1.name()+"&"+f2.name()+")",
            f1.fvMsh()
        )
    );

    staggeredScalarField& Prod = tProd.ref();

    Prod = Zero;

    forAllCells(Prod, d1, i, j, k)
        for (int d2 = 0; d2 < 3; d2++)
            Prod(d1,i,j,k) +=
                f1(d1,i,j,k)[d2]
              * (
                    d1 == d2
                  ? f2(d2,i,j,k)[d1]
                  : 0.25
                  * (
                        f2(d2,i,j,k)[d1]
                      + f2(d2,lowerNei(i,j,k,d1))[d1]
                      + f2(d2,upperNei(i,j,k,d2))[d1]
                      + f2(d2,lowerNei(upperNei(i,j,k,d2),d1))[d1]
                    )
                );

    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const tmp<staggeredVectorField>& tf1,
    const staggeredVectorField& f2
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(tf1(),f2));
    tf1.clear();
    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const staggeredVectorField& f1,
    const tmp<staggeredVectorField>& tf2
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(f1,tf2()));
    tf2.clear();
    return tProd;
}

inline tmp<staggeredScalarField>
stagDotProduct
(
    const tmp<staggeredVectorField>& tf1,
    const tmp<staggeredVectorField>& tf2
)
{
    tmp<staggeredScalarField> tProd(stagDotProduct(tf1(),tf2()));
    tf1.clear();
    tf2.clear();
    return tProd;
}

}

}

}

#endif
