#ifndef StencilSpace_H
#define StencilSpace_H

#include "CellSpace.H"
#include "scalar.H"

namespace Foam
{

template<class Form, direction Ncmpts>
class StencilSpace
:
    public CellSpace<Form, scalar, Ncmpts>
{

public:

    // Stencil space type
    typedef StencilSpace<Form, Ncmpts> ssType;

    // Cell space type
    typedef CellSpace<Form, scalar, Ncmpts> csType;

    // Component type
    typedef scalar cmpt;

    // Component of component type
    typedef scalar cmptCmpt;

    // Construct null
    inline StencilSpace();

    // Construct initialised to zero
    inline StencilSpace(const Foam::zero);

    // Construct uniform
    inline StencilSpace(const scalar&);

    // Construct from cell space
    inline StencilSpace(const CellSpace<Form, scalar, Ncmpts>&);

    // Construct from Istream
    inline StencilSpace(Istream&);
};

// Stencil space should always reduce to self

template<class Form, class Cmpt, direction Ncmpts>
class CellSpaceCast<StencilSpace<Form,Ncmpts>, StencilSpace<Form,Ncmpts>, Cmpt>
{
public:
    typedef Form type;
};

// Class to cast two stencil forms into one. The highest order stencil should
// always prevail. Must be specialized by each combination of stencils.

template<class Type1, class Type2>
class StencilSpaceCast
{};

// By default type T is not a stencil. Stencils should specialize.

template<class T>
inline bool isStencil()
{
    return false;
}

}

#include "StencilSpaceI.H"

#endif