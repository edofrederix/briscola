#ifndef linearProlongationScheme_H
#define linearProlongationScheme_H

#include "prolongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearProlongationScheme
:
    public prolongationScheme<Type,MeshType>
{
    // Cached weights
    meshField<tensor,MeshType> weights_;

    // Set weights
    void setWeights();

protected:

    // MidPoint prolongation function
    template<template<class> class OpType>
    void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const OpType<Type>&
    );

public:

    // Runtime type information
    TypeName("linear");

    // Construct from dictionary and mesh
    linearProlongationScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearProlongationScheme(const fvMesh&);

    // In-place prolongation from a coarse cell level to the adjacent fine cell level,
    // using the = operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const eqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell level,
    // using the += operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const plusEqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell level,
    // using the -= operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const minusEqOp<Type>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearProlongationScheme.C"
#endif

#endif
