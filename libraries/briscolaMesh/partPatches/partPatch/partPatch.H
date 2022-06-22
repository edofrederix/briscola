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
    // Reference to the mesh
    const mesh& msh_;

    // Copy of the dictionary
    dictionary dict_;

    // Name of the patch
    const word name_;

    // Offset of this patch to the boundary
    labelVector boundaryOffset_;

    // Rotation of the neighbor in case of a parallel patch (defaults to the
    // identity tensor)
    labelTensor T_;

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

    // Selector
    static autoPtr<partPatch> New
    (
        const mesh& msh,
        const dictionary& dict
    );

    // Destructor
    virtual ~partPatch();

    // Access to the mesh
    inline const mesh& msh() const { return msh_; }

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

    // Rotation of the neighbor
    inline const labelTensor& T() const { return T_; }
};

}

}

#endif
