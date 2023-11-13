#ifndef partPatch_H
#define partPatch_H

#include "dictionary.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "Pstream.H"

namespace Foam
{

namespace briscola
{

class mesh;

class partPatch
{
protected:

    // Copy of the dictionary
    dictionary dict_;

    // Name of the patch
    const word name_;

    // Boundary offset of this patch
    labelVector boundaryOffset_;

    // Boundary offset degree, giving the number on non-zero components in the
    // boundary offset vector
    label boundaryOffsetDegree_;

    // Rotation of the neighbor in case of a parallel patch (defaults to the
    // identity tensor)
    labelTensor T_;

    // Master bool. Part patches that require communication can be master or
    // slave. By default a part patch is master. Can be modified by derived
    // classes.
    bool master_;

public:

    // Runtime type information
    TypeName("partPatch");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        partPatch,
        dictionary,
        (
            const mesh& msh,
            const dictionary& dict
        ),
        (msh, dict)
    );

    // Construct from components
    partPatch(const mesh& msh, const dictionary& dict);

    // Copy constructor
    partPatch(const partPatch& pp);

    // Clone part patch
    virtual autoPtr<partPatch> clone() const = 0;

    // Selector
    static autoPtr<partPatch> New
    (
        const mesh& msh,
        const dictionary& dict
    );

    // Destructor
    virtual ~partPatch();

    // Access to the dictionary
    inline dictionary& dict() { return dict_; }
    inline const dictionary& dict() const { return dict_; }

    // Name of the patch
    inline const word& name() const { return name_; }

    // Offset of the patch to the boundary
    inline const labelVector& boundaryOffset() const
    {
        return boundaryOffset_;
    }

    // Boundary offset degree
    inline label boundaryOffsetDegree() const
    {
        return boundaryOffsetDegree_;
    }

    // Rotation of the neighbor
    inline const labelTensor& T() const { return T_; }

    // Master bool
    inline bool master() const { return master_; }
    inline bool slave() const { return !master_; }

    // Part patch type number
    virtual label typeNum() const = 0;

    // Neighbor processor number
    virtual label neighborProcNum() const { return -1; }
};

}

}

#endif
