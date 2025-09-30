#ifndef boundary_H
#define boundary_H

#include "dictionary.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "Pstream.H"
#include "faceLabel.H"

namespace Foam
{

namespace briscola
{

class mesh;

class boundary
{
protected:

    // Copy of the dictionary
    dictionary dict_;

    // Name of the boundary
    const word name_;

    // Boundary offset of this boundary
    labelVector offset_;

    // Boundary offset degree, giving the number on non-zero components in the
    // boundary offset vector
    label offsetDegree_;

    // Master bool. Boundaries that require communication can be master or
    // slave. By default a boundary is master. Can be modified by derived
    // classes.
    bool master_;

    // Extension of the boundary into other boundaries
    faceLabel extension_;

    // Function to make the boundary an extended one
    void extend();

public:

    // Runtime type information
    TypeName("boundary");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        boundary,
        dictionary,
        (
            const mesh& msh,
            const dictionary& dict
        ),
        (msh, dict)
    );

    // Construct from components
    boundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    boundary(const boundary&);

    // Clone part boundary
    virtual autoPtr<boundary> clone() const = 0;

    // Selector
    static autoPtr<boundary> New
    (
        const mesh& msh,
        const dictionary& dict
    );

    // Destructor
    virtual ~boundary();

    // Check if cast to derived type is possible
    template<class BoundaryType>
    inline bool castable() const
    {
        return dynamic_cast<const BoundaryType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class BoundaryType>
    inline BoundaryType& cast()
    {
        BoundaryType* bPtr = dynamic_cast<BoundaryType*>(this);

        if (!bPtr)
        {
            FatalErrorInFunction
                << "Boundary " << name_ << " cannot be cast to "
                << BoundaryType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *bPtr;
    }

    // Cast to derived const type, if possible. Returns error otherwise.
    template<class BoundaryType>
    inline const BoundaryType& cast() const
    {
        return const_cast<boundary*>(this)->cast<BoundaryType>();
    }

    // Access to the dictionary
    inline dictionary& dict() { return dict_; }
    inline const dictionary& dict() const { return dict_; }

    // Name of the boundary
    inline const word& name() const { return name_; }

    // Offset of the boundary
    inline const labelVector& offset() const
    {
        return offset_;
    }

    // Boundary offset degree
    inline label offsetDegree() const
    {
        return offsetDegree_;
    }

    // Master bool
    inline bool master() const { return master_; }
    inline bool slave() const { return !master_; }

    // Access to the extension face label
    inline faceLabel& extension() { return extension_; }
    inline const faceLabel& extension() const { return extension_; }

    // Check if the boundary is extended
    inline bool extended() const { return extension_ != faceLabel::zero; }

    // Neighbor processor number
    virtual label neighborProcNum() const { return -1; }
};

}

}

#endif
