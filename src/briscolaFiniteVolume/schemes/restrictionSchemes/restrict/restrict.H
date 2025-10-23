#ifndef restrict_H
#define restrict_H

#include "meshFields.H"
#include "faceFields.H"
#include "colocated.H"
#include "staggered.H"

#include "restrictionScheme.H"
#include "faceRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Restriction

template<class Type, class MeshType>
void restrict(const meshField<Type,MeshType>& field)
{
    const_cast<meshField<Type,MeshType>&>(field).makeDeep();

    return restrictionScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        field.restrictionScheme()
    )->restrict(const_cast<meshField<Type,MeshType>&>(field));
}

template<class Type, class MeshType>
void restrict(const faceField<Type,MeshType>& field)
{
    forAll(field, fd)
        const_cast<meshField<Type,MeshType>&>(field[fd]).makeDeep();

    return faceRestrictionScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        field.restrictionScheme()
    )->restrict(const_cast<faceField<Type,MeshType>&>(field));
}

// Collapse fields by making them shallow

template<class Type, class MeshType>
void collapse(const meshField<Type,MeshType>& field)
{
    const_cast<meshField<Type,MeshType>&>(field).makeShallow();
}

template<class Type, class MeshType>
void collapse(const faceField<Type,MeshType>& field)
{
    forAll(field, s)
        collapse(field[s]);
}

}

}

}

#endif
