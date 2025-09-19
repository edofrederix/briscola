#ifndef patch_H
#define patch_H

#include "meshObject.H"
#include "face.H"
#include "NamedEnum.H"

namespace Foam
{

namespace briscola
{

class geometry;

class patch
:
    public meshObject<geometry>
{
public:

    // Enumeration listing the patch types
    enum patchType
    {
        PATCH,
        PERIODIC,
        MAPPED,
        EMPTY
    };

    static const NamedEnum<patchType,4> patchTypeNames;

private:

    // Copy of the patch dictionary
    const dictionary dict_;

    // Name of the patch
    const word name_;

    // Patch type
    const patchType type_;

    // List of pointers to the faces contained in this patch
    List<const face*> facePtrs_;

    // Check the consistency of this patch
    void checkConsistency() const;

public:

    // Construct from dictionary
    patch(const geometry&, const label, const word, const dictionary&);

    // Construct from list of pointers to faces
    patch(const geometry&, const label, const word, const List<const face*>&);

    // Copy constructor
    patch(const patch&);

    // Clone this patch (needed by PtrList)
    inline autoPtr<patch> clone() const
    {
        return autoPtr<patch>(new patch(*this));
    }

    // Destructor
    virtual ~patch();

    // Access to the patch name
    inline word name() const
    {
        return name_;
    }

    patchType type() const
    {
        return type_;
    }

    word typeName() const
    {
        return patchTypeNames[type_];
    }

    // Access to the parent geometry
    inline const geometry& parentGeometry() const
    {
        return parent_;
    }

    // Access to the patch dictionary
    inline const dictionary& dict() const
    {
        return dict_;
    }

    // Access to the list of face pointers
    inline const List<const face*>& facePtrs() const
    {
        return facePtrs_;
    }

    // Check if given face is part of this patch
    inline bool hasFace(const face& f) const
    {
        forAll(facePtrs_, i)
        {
            if (*facePtrs_[i] == f) return true;
        }

        return false;
    }

    // Equality operator
    bool operator==(const patch&) const;

    // Unequality operator
    bool operator!=(const patch&) const;
};

}

}

#endif
