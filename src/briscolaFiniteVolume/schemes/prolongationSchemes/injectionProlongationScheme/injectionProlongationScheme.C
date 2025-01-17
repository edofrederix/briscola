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
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const labelVector R(coarse.mshPart().R());

    forAllCells(fine, i, j, k)
    {
        const label il = i/R.x();
        const label jl = j/R.y();
        const label kl = k/R.z();

        bop
        (
            fine(i,j,k),
            coarse(il,jl,kl)
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
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void injectionProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void injectionProlongationScheme<Type,MeshType>::prolong
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
