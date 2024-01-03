#ifndef exSchemesReconstruct_H
#define exSchemesReconstruct_H

#include "meshField.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Reconstruct a colocated field from a staggered one

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<Type,staggered>& field
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

    const meshField<faceVector,colocated>& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    forAllCells(Recon, i, j, k)
        Recon(i,j,k) =
            0.5
          * (
              - field(0,i,  j,  k  ) * fn(i,j,k).left()
              + field(0,i+1,j,  k  ) * fn(i,j,k).right()
              - field(1,i,  j,  k  ) * fn(i,j,k).bottom()
              + field(1,i,  j+1,k  ) * fn(i,j,k).top()
              - field(2,i,  j,  k  ) * fn(i,j,k).aft()
              + field(2,i,  j,  k+1) * fn(i,j,k).fore()
            );

    return tRecon;
}

// inline tmp<colocatedVectorField>
// reconstruct
// (
//     const colocatedLowerFaceScalarField& field
// )
// {

//     tmp<colocatedVectorField> tRecon
//     (
//         new colocatedVectorField
//         (
//             "reconstruct("+field.name()+")",
//             field.fvMsh()
//         )
//     );

//     colocatedVectorField& Recon = tRecon.ref();

//     const meshField<faceVector,colocated>& fn =
//         field.fvMsh().template metrics<colocated>().faceNormals();

//     const meshField<faceScalar,colocated>& fa =
//         field.fvMsh().template metrics<colocated>().faceAreas();

//     forAllCells(Recon, i, j, k)
//     {
//         vector totalWeight =
//                   - fa(i,j,k).left()   * fn(i,j,k).left()
//                   + fa(i,j,k).right()  * fn(i,j,k).right()
//                   - fa(i,j,k).bottom() * fn(i,j,k).bottom()
//                   + fa(i,j,k).top()    * fn(i,j,k).top()
//                   - fa(i,j,k).aft()    * fn(i,j,k).aft()
//                   + fa(i,j,k).fore()   * fn(i,j,k).fore();

//         vector totalSum =
//                   - field(i,  j,  k  ).left()   * fn(i,j,k).left()
//                   + field(i+1,j,  k  ).left()   * fn(i,j,k).right()
//                   - field(i,  j,  k  ).bottom() * fn(i,j,k).bottom()
//                   + field(i,  j+1,k  ).bottom() * fn(i,j,k).top()
//                   - field(i,  j,  k  ).aft()    * fn(i,j,k).aft()
//                   + field(i,  j,  k+1).aft()    * fn(i,j,k).fore();

//         Recon(i,j,k) =
//             vector
//             (
//                 totalSum[0] / totalWeight[0],
//                 totalSum[1] / totalWeight[1],
//                 totalSum[2] / totalWeight[2]
//             );
//     }

//     return tRecon;
// }

inline tmp<staggeredScalarField>
transposeMultiplicate
(
    const staggeredVectorField& field1,
    const staggeredVectorField& field2
)
{

    tmp<staggeredScalarField> tRecon
    (
        new staggeredScalarField
        (
            "("+field1.name()+")",
            field1.fvMsh()
        )
    );

    staggeredScalarField& Recon = tRecon.ref();

    forAllCells(Recon, d, i, j, k)
    {
        Recon(d,i,j,k) =
                field1(d,i,j,k)[0] * field2(0,i,j,k)[d]
              + field1(d,i,j,k)[1] * field2(1,i,j,k)[d]
              + field1(d,i,j,k)[2] * field2(2,i,j,k)[d];
    }

    return tRecon;
}

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<LowerFaceSpace<Type>,colocated>& field
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

    const meshField<faceVector,colocated>& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    const meshField<faceScalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(Recon, i, j, k)
    {
        vector totalWeight =
                  - fa(i,j,k).left()   * fn(i,j,k).left()
                  + fa(i,j,k).right()  * fn(i,j,k).right()
                  - fa(i,j,k).bottom() * fn(i,j,k).bottom()
                  + fa(i,j,k).top()    * fn(i,j,k).top()
                  - fa(i,j,k).aft()    * fn(i,j,k).aft()
                  + fa(i,j,k).fore()   * fn(i,j,k).fore();

        ReconType totalSum =
                + field(i,  j,  k  ).left()   * fn(i,j,k).left()
                - field(i+1,j,  k  ).left()   * fn(i,j,k).right()
                + field(i,  j,  k  ).bottom() * fn(i,j,k).bottom()
                - field(i,  j+1,k  ).bottom() * fn(i,j,k).top()
                + field(i,  j,  k  ).aft()    * fn(i,j,k).aft()
                - field(i,  j,  k+1).aft()    * fn(i,j,k).fore();

        Recon(i,j,k).x() = totalSum.x() / totalWeight.x();
        Recon(i,j,k).y() = totalSum.y() / totalWeight.y();
        Recon(i,j,k).z() = totalSum.z() / totalWeight.z();
    }

    return tRecon;
}

}

}

}

}

#endif
