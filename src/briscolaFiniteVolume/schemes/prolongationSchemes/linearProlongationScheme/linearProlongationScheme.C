#include "linearProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::setWeightsAndIndices()
{
    weights_.clear();
    indices_.clear();

    weights_.setSize(this->fvMsh().size()-1);
    indices_.setSize(this->fvMsh().size()-1);

    for (label l = 0; l < this->fvMsh().size()-1; l++)
    {
        const labelVector R(this->fvMsh()[l+1].R());
        const labelVector offset(this->fvMsh()[l+1].aggParentOffset());

        weights_[l].setSize(MeshType::numberOfDirections);
        indices_[l].setSize(MeshType::numberOfDirections);

        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            const meshDirection<vector,MeshType>& ccf =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l][d];

            const meshDirection<vector,MeshType>& ccc =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l+1][d];

            List<FixedList<scalar,8>>& weights = weights_[l][d];
            List<FixedList<labelVector,8>>& indices = indices_[l][d];

            weights.setSize(ccf.size());
            indices.setSize(ccf.size());

            label c = 0;
            forAllCells(ccf, i, j, k)
            {
                const labelVector ijk(i,j,k);
                const labelVector off
                (
                    i%2 ? 1 : -1,
                    j%2 ? 1 : -1,
                    k%2 ? 1 : -1
                );

                // Set the indices accounting for possible aggregate parent
                // offset

                int v = 0;
                labelVector abc;
                for (abc.z() = 0; abc.z() < 2; abc.z()++)
                    for (abc.y() = 0; abc.y() < 2; abc.y()++)
                        for (abc.x() = 0; abc.x() < 2; abc.x()++)
                            indices[c][v++] =
                                briscola::cmptDivide(ijk,R)
                              + briscola::cmptMultiply(abc,off)
                              + offset;

                // Interpolation cell vertices

                vertexVector vertices;
                for (v = 0; v < 8; v++)
                    vertices[v] = ccc(indices[c][v]);

                // Calculate the tri-linear interpolation point using a large
                // tolerance. As tolerance length scale the cubic root of the
                // interpolation hexahedron volume is used.

                vector w =
                    interpolationWeights
                    (
                        ccf(i,j,k),
                        vertices,
                        true,
                        false,
                        1e-3*Foam::cbrt(Foam::mag(hexVolume(vertices)))
                    );

                // On very coarse unstructured meshes the interpolation point
                // calculation may fail. In that case we take the midpoint.

                if (w == -vector::one)
                    w = vector::one/2.0;

                // Set the weights

                v = 0;
                for (abc.z() = 0; abc.z() < 2; abc.z()++)
                    for (abc.y() = 0; abc.y() < 2; abc.y()++)
                        for (abc.x() = 0; abc.x() < 2; abc.x()++)
                            weights[c][v++] =
                                (1.0 - abc.x() + w.x()*(2*abc.x() - 1))
                              * (1.0 - abc.y() + w.y()*(2*abc.y() - 1))
                              * (1.0 - abc.z() + w.z()*(2*abc.z() - 1));

                c++;
            }
        }
    }
}

template<class Type, class MeshType>
template<template<class> class OpType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const label l = fine.levelNum();

    const List<List<FixedList<scalar,8>>>& weights = weights_[l];
    const List<List<FixedList<labelVector,8>>>& indices = indices_[l];

    // Linear prolongation uses ghost cells, so boundaries need to be corrected

    const_cast<meshLevel<Type,MeshType>&>(coarse).correctBoundaryConditions();
    const_cast<meshLevel<Type,MeshType>&>(coarse).correctAggData();

    forAll(fine, d)
    {
        label c = 0;
        forAllCells(fine[d], i, j, k)
        {
            Type value = Zero;

            for (label v = 0; v < 8; v++)
                if (weights[d][c][v] != 0.0)
                    value += weights[d][c][v]*coarse(d,indices[d][c][v]);

            bop(fine(d,i,j,k), value);

            c++;
        }
    }
}

template<class Type, class MeshType>
linearProlongationScheme<Type,MeshType>::linearProlongationScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    prolongationScheme<Type,MeshType>(fvMsh, is)
{
    setWeightsAndIndices();
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const minusEqOp<Type>& bop
)
{
    this->prolong<minusEqOp>(fine, coarse, bop);
}

}

}

}
