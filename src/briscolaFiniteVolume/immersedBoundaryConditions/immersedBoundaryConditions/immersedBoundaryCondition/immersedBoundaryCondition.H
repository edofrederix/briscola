#ifndef immersedBoundaryCondition_H
#define immersedBoundaryCondition_H

#include "meshFieldFwd.H"
#include "immersedBoundaryConditionFwd.H"
#include "immersedBoundaryFwd.H"
#include "rectilinearMesh.H"

#include "List.H"
#include "Switch.H"
#include "dictionary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class Type, class MeshType>
class immersedBoundaryCondition
:
    public refCount
{
protected:

    // Reference to mesh
    const fvMesh& fvMsh_;

    // Reference to the mesh field
    meshField<Type,MeshType>& field_;

    // Reference to immersed boundary
    const immersedBoundary<MeshType>& ib_;

    // Pointer to the forcing mask
    const meshField<label,MeshType>* forcingMaskPtr_;

    // Copy of the dictionary
    dictionary dict_;

    // Relaxation factor
    scalar omega_;

    // Shape overlap warning
    static void shapeOverlapWarning()
    {
        WarningInFunction
            << "Overlapping shapes identified. This may cause issues "
            << "with " << typeName << " immersed boundary condition."
            << endl;
    }

    // Function to read boundary value entry from the dictionary
    List<Type> read(const word) const;

public:

    // Declare name of the class and its debug switch
    TypeName("immersedBoundaryCondition");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        immersedBoundaryCondition,
        dictionary,
        (
            const meshField<Type,MeshType>& field,
            const immersedBoundary<MeshType>& ib
        ),
        (field, ib)
    );

    // Construct from field, immersed boundary and mask pointer
    immersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib,
        const meshField<label,MeshType>* maskPtr
    );

    // Copy constructor
    immersedBoundaryCondition
    (
        const immersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    immersedBoundaryCondition
    (
        const immersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Selector from dictionary
    static autoPtr<immersedBoundaryCondition<Type,MeshType>> New
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~immersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const = 0;

    // Construct and return a clone for alternative field
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>&
    ) const = 0;

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label);
    virtual void evaluate(const label, const label) = 0;

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the mesh field
    inline const meshField<Type,MeshType>& field() const
    {
        return field_;
    }

    // Access to the immersed boundary
    inline const immersedBoundary<MeshType>& ib() const
    {
        return ib_;
    }

    // Const access to the forcing mask pointer
    inline const meshField<label,MeshType>* forcingMaskPtr() const
    {
        return forcingMaskPtr_;
    }

    // Const access to the forcing mask (unchecked!)
    inline const meshField<label,MeshType>& forcingMask() const
    {
        return *forcingMaskPtr_;
    }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }
};

}

}

}

#ifdef NoRepository
    #include "immersedBoundaryCondition.C"
#endif

#endif
