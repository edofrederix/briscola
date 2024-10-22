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

// Staggered transpose, exchanging vector directions with mesh directions

inline tmp<staggeredVectorField> stagT(const staggeredVectorField& f)
{
    const_cast<staggeredVectorField&>(f).correctBoundaryConditions();

    tmp<staggeredVectorField> tT
    (
        new staggeredVectorField
        (
            "T("+f.name()+")",
            f.fvMsh()
        )
    );

    staggeredVectorField& T = tT.ref();

    T = Zero;

    forAllCells(T, di, i, j, k)
        for (int dj = 0; dj < 3; dj++)
            T(di,i,j,k)[dj] =
                di == dj
              ? f(dj,i,j,k)[di]
              : 0.25
              * (
                    f(dj,i,j,k)[di]
                  + f(dj,lowerNei(i,j,k,di))[di]
                  + f(dj,upperNei(i,j,k,dj))[di]
                  + f(dj,lowerNei(upperNei(i,j,k,dj),di))[di]
                );

    return tT;
}

inline tmp<staggeredVectorField> stagT(const tmp<staggeredVectorField>& tf)
{
    tmp<staggeredVectorField> tT(stagT(tf()));
    tf.clear();
    return tT;
}

}

}

}

#endif
