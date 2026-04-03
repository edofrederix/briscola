#include "injectionProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
template<template<class> class OpType>
void injectionProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const labelVector R(coarse.lvl().R());
    const labelVector offset(coarse.lvl().aggParentOffset());

    const_cast<meshLevel<Type,MeshType>&>(coarse).correctAggData();

    forAllCells(fine, d, i, j, k)
    {
        bop
        (
            fine(d,i,j,k),
            coarse(d, briscola::cmptDivide(labelVector(i,j,k),R) + offset)
        );
    }
}

template<class Type, class MeshType>
injectionProlongationScheme<Type,MeshType>::injectionProlongationScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    prolongationScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
void injectionProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void injectionProlongationScheme<Type,MeshType>::prolong
(
    meshLevel<Type,MeshType>& fine,
    const meshLevel<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void injectionProlongationScheme<Type,MeshType>::prolong
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
