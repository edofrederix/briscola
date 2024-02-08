#ifndef immersedBoundaryCondition_H
#define immersedBoundaryCondition_H

#include "immersedBoundaryConditionFwd.H"
#include "immersedBoundaryFwd.H"

#include "List.H"
#include "shape.H"
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

    // Copy of the dictionary
    dictionary dict_;

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

    // Construct from solverDict and fvMesh
    immersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Selector from dictionary
    static autoPtr<immersedBoundaryCondition<Type,MeshType>> New
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~immersedBoundaryCondition();

    // Correct the stencils of a linear system to apply the IBM correction.
    // This function is overwritten by fully implicit IBM types and does
    // nothing otherwise.
    virtual void correctLinearSystem
    (
        linearSystem<stencil,Type,MeshType>&
    );

    // Correct the forcing points as Jacobi points. This function is
    // overwritten by IBM types with large stencils (> 7 points in 3D).
    virtual void correctJacobiPoints
    (
        meshLevel<Type,MeshType>&
    );

    // Return colocated volume source term for a staggered velocity field.
    // Only defined for staggeredScalarImmersedBoundary
    // tmp<meshField<scalar,colocated>> IBMSource
    // (
    //     const meshField<scalar,staggered>&
    // );

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
