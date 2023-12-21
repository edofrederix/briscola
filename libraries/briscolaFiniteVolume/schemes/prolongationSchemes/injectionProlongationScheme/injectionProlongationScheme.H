#ifndef injectionProlongationScheme_H
#define injectionProlongationScheme_H

#include "prolongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class injectionProlongationScheme
:
    public prolongationScheme<Type,MeshType>
{
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
    TypeName("injection");

    // Construct from dictionary and mesh
    injectionProlongationScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    injectionProlongationScheme(const fvMesh&);

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the = operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const eqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the += operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const plusEqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the -= operation
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
    #include "injectionProlongationScheme.C"
#endif

#endif
