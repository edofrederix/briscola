#ifndef scheme_H
#define scheme_H

#include "fvMesh.H"
#include "dictionary.H"
#include "ops.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Function to check of OpType is eqOp
template<template<class> class OpType>
inline bool isEqOp() { return false; }

// Specialization
template<>
inline bool isEqOp<eqOp>() { return true; }

// Function to check of OpType is eqMinusOp
template<template<class> class OpType>
inline bool isEqMinusOp() { return false; }

// Specialization
template<>
inline bool isEqMinusOp<eqMinusOp>() { return true; }

// Flip operator

template<class Type, template<class> class OpType>
class FlipOp
{
public:
    typedef OpType<Type> type;
};

template<class Type>
class FlipOp<Type,eqOp>
{
public:
    typedef eqMinusOp<Type> type;
};

template<class Type>
class FlipOp<Type,eqMinusOp>
{
public:
    typedef eqOp<Type> type;
};

template<class Type>
class FlipOp<Type,plusEqOp>
{
public:
    typedef minusEqOp<Type> type;
};

template<class Type>
class FlipOp<Type,minusEqOp>
{
public:
    typedef plusEqOp<Type> type;
};

class scheme
{
    // Copy of the scheme dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

public:

    // Runtime type information
    TypeName("scheme");

    // Construct from dictionary
    scheme(const dictionary&, const fvMesh& fvMsh);

    // Copy constructor
    scheme(const scheme&);

    // Destructor
    virtual ~scheme();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }
};

}

}

}

#endif
