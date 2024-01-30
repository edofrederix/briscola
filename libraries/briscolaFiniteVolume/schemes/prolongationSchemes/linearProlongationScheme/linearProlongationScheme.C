#include "linearProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::setWeights()
{
    for (label l = 0; l < this->fvMsh().size()-1; l++)
    {
        const labelVector R(this->fvMsh()[l+1].R());

        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            const meshDirection<vector,MeshType>& ccf =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l][d];

            const meshDirection<vector,MeshType>& ccc =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l+1][d];

            meshDirection<vertexScalar,MeshType>& weights = weights_[l][d];

            weights = Zero;

            forAllCells(ccf, i, j, k)
            {
                const label ox = (i%2)*2-1;
                const label oy = (j%2)*2-1;
                const label oz = (k%2)*2-1;

                vertexVector vertices;

                int q = 0;
                for (int c = 0; c < 2; c++)
                    for (int b = 0; b < 2; b++)
                        for (int a = 0; a < 2; a++)
                            vertices[q++] =
                                ccc
                                (
                                    i/R.x() + a*ox,
                                    j/R.y() + b*oy,
                                    k/R.z() + c*oz
                                );

                // Calculate the tri-linear interpolation point. On very coarse
                // unstructured meshes this may fail. In that case we take the
                // mid point.

                vector w =
                    interpolationWeights(ccf(i,j,k), vertices, true, false);

                if (w == -vector::one)
                    w = vector::one/2.0;

                // Set the weights

                weights(i,j,k) =
                    vertexScalar
                    (
                        (1.0-w.x()) * (1.0-w.y()) * (1.0-w.z()),
                        (    w.x()) * (1.0-w.y()) * (1.0-w.z()),
                        (1.0-w.x()) * (    w.y()) * (1.0-w.z()),
                        (    w.x()) * (    w.y()) * (1.0-w.z()),
                        (1.0-w.x()) * (1.0-w.y()) * (    w.z()),
                        (    w.x()) * (1.0-w.y()) * (    w.z()),
                        (1.0-w.x()) * (    w.y()) * (    w.z()),
                        (    w.x()) * (    w.y()) * (    w.z())
                    );
            }
        }
    }
}

template<class Type, class MeshType>
template<template<class> class OpType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const labelVector R(coarse.mshPart().R());

    meshDirection<vertexScalar,MeshType>& weights =
        weights_[fine.levelNum()][fine.directionNum()];

    forAllCells(fine, i, j, k)
    {
        const label ox = ((i%2)*2-1) * (R.x()>1);
        const label oy = ((j%2)*2-1) * (R.y()>1);
        const label oz = ((k%2)*2-1) * (R.z()>1);

        Type value = Zero;

        int q = 0;
        for (int c = 0; c < 2; c++)
            for (int b = 0; b < 2; b++)
                for (int a = 0; a < 2; a++)
                    value +=
                        weights(i,j,k)[q++]
                      * coarse
                        (
                            i/R.x() + a*ox,
                            j/R.y() + b*oy,
                            k/R.z() + c*oz
                        );

        bop(fine(i,j,k), value);
    }
}

template<class Type, class MeshType>
linearProlongationScheme<Type,MeshType>::linearProlongationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dict,fvMsh),
    weights_
    (
        "weights",
        this->fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    )
{
    setWeights();
}

template<class Type, class MeshType>
linearProlongationScheme<Type,MeshType>::linearProlongationScheme
(
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dictionary(),fvMsh),
    weights_
    (
        "weights",
        this->fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    )
{
    setWeights();
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const minusEqOp<Type>& bop
)
{
    this->prolong<minusEqOp>(fine, coarse, bop);
}

}

}

}
