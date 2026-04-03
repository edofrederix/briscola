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
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const OpType<Type>&
    );

public:

    // Runtime type information
    TypeName("injection");

    // Construct from stream and mesh
    injectionProlongationScheme(const fvMesh&, Istream& = scheme::nullStream);

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the = operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const eqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the += operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const plusEqOp<Type>&
    );

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the -= operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
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
