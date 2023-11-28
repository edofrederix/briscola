#include "midPointProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
template<template<class> class OpType>
void midPointProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const labelVector R(coarse.mshPart().R());
    const vector s(MeshType::shift[coarse.directionNum()]);

    forAllCells(fine, i, j, k)
    {
        const label il = i/R.x();
        const label jl = j/R.y();
        const label kl = k/R.z();

        const label ox = (i % 2)*2 - 1;
        const label oy = (j % 2)*2 - 1;
        const label oz = (k % 2)*2 - 1;

        const label iu = il + (R.x() == 1 ? 0 : ox);
        const label ju = jl + (R.y() == 1 ? 0 : oy);
        const label ku = kl + (R.z() == 1 ? 0 : oz);

        bop
        (
            fine(i,j,k),
            (
                (1.5+ox*s.x())*(1.5+oy*s.y())*(1.5+oz*s.z())*coarse(il,jl,kl)
              + (0.5-ox*s.x())*(1.5+oy*s.y())*(1.5+oz*s.z())*coarse(iu,jl,kl)
              + (1.5+ox*s.x())*(0.5-oy*s.y())*(1.5+oz*s.z())*coarse(il,ju,kl)
              + (0.5-ox*s.x())*(0.5-oy*s.y())*(1.5+oz*s.z())*coarse(iu,ju,kl)
              + (1.5+ox*s.x())*(1.5+oy*s.y())*(0.5-oz*s.z())*coarse(il,jl,ku)
              + (0.5-ox*s.x())*(1.5+oy*s.y())*(0.5-oz*s.z())*coarse(iu,jl,ku)
              + (1.5+ox*s.x())*(0.5-oy*s.y())*(0.5-oz*s.z())*coarse(il,ju,ku)
              + (0.5-ox*s.x())*(0.5-oy*s.y())*(0.5-oz*s.z())*coarse(iu,ju,ku)
            )
          / 8.0
        );
    }
}

template<class Type, class MeshType>
midPointProlongationScheme<Type,MeshType>::midPointProlongationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointProlongationScheme<Type,MeshType>::midPointProlongationScheme
(
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
void midPointProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void midPointProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void midPointProlongationScheme<Type,MeshType>::prolong
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
