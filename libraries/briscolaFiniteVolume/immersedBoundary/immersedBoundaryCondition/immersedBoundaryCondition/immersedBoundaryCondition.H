#ifndef immersedBoundaryCondition_H
#define immersedBoundaryCondition_H

#include "immersedBoundaryConditionFwd.H"
#include "immersedBoundaryFwd.H"

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
class meshField;

template<class SType, class Type, class MeshType>
class linearSystem;

template<class Type, class MeshType>
class immersedBoundaryCondition
{
protected:

    // Reference to mesh
    const fvMesh& fvMsh_;

    // Reference to immersed boundary
    const immersedBoundary<MeshType>& IB_;

    // Reference to the forcing points mask
    const meshField<label,MeshType>& forcingPoints_;

    // Copy of the dictionary
    dictionary dict_;

    // Relaxation factor
    scalar omega_;

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
            const meshField<Type,MeshType>& mshField,
            const immersedBoundary<MeshType>& ib
        ),
        (mshField, ib)
    );

    // Construct from meshField and immersedBoundary
    immersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib,
        const meshField<label,MeshType>& fp
    );

    // Selector from dictionary
    static autoPtr<immersedBoundaryCondition<Type,MeshType>> New
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~immersedBoundaryCondition();

    // Correct the forcing points as Jacobi points.
    virtual void correctJacobiPoints
    (
        meshLevel<Type,MeshType>&
    ) const = 0;

    // Inline access functions
    inline const dictionary& dict() const { return dict_; }

    inline const meshField<label,MeshType>& forcingPoints() const
    { return forcingPoints_; }

};

}

}

}

#ifdef NoRepository
    #include "immersedBoundaryCondition.C"
#endif

#endif
