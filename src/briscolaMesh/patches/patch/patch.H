#ifndef patch_H
#define patch_H

#include "meshObject.H"
#include "face.H"

namespace Foam
{

namespace briscola
{

class geometry;

class patch
:
    public meshObject<geometry>
{
private:

    // Copy of the patch dictionary
    const dictionary dict_;

    // Name of the patch
    const word name_;

    // List of pointers to the faces contained in this patch
    List<const face*> facePtrs_;

    // Check the consistency of this patch
    void checkConsistency() const;

public:

    // Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        patch,
        dictionary,
        (
            const geometry& geo,
            const label num,
            const word name,
            const dictionary& dict
        ),
        (geo, num, name, dict)
    );

    // Construct from dictionary
    patch(const geometry&, const label, const word, const dictionary&);

    // Copy constructor
    patch(const patch&);

    // Clone this patch
    inline autoPtr<patch> clone() const
    {
        return autoPtr<patch>(new patch(*this));
    }

    // Selector
    static autoPtr<patch> New
    (
        const geometry& geo,
        const label num,
        const word name,
        const dictionary& dict
    );

    // Destructor
    virtual ~patch();

    // Check if cast to derived type is possible
    template<class PatchType>
    inline bool castable() const
    {
        return dynamic_cast<const PatchType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class PatchType>
    inline PatchType& cast()
    {
        PatchType* ptr = dynamic_cast<PatchType*>(this);

        if (!ptr)
        {
            FatalErrorInFunction
                << "Patch cannot be cast to " << PatchType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *ptr;
    }

    // Cast to derived const type, if possible. Returns error otherwise.
    template<class PatchType>
    inline const PatchType& cast() const
    {
        const PatchType* ptr = dynamic_cast<const PatchType*>(this);

        if (!ptr)
        {
            FatalErrorInFunction
                << "Patch cannot be cast to " << PatchType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *ptr;
    }

    // Access to the patch name
    inline word name() const
    {
        return name_;
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
