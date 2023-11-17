#ifndef immersedBoundaryMethod_H
#define immersedBoundaryMethod_H

#include "immersedBoundaryMethodFwd.H"

#include "List.H"
#include "shape.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvMeshMetrics.H"
#include "meshField.H"
#include "linearSystem.H"
#include "rectilinearMesh.H"
#include "briscolaUtils.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class immersedBoundaryMethod
{
protected:
    // Reference to mesh
    const fvMesh& fvMsh_;

    // Mass source active switch
    bool massSourceActive_;

    // True for ghost-cell methods updating values as Jacobi points
    bool JacobiGhostMethod_;

    // Convergence relaxation parameter
    scalar omega_;

    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<label,MeshType> mask_;

    // The ghost cell mask is equal to 1 in IB cells with at least
    // one neighboring fluid cell and 0 in other cells
    meshField<label,MeshType> ghostMask_;

    // List of shapes making up the immersed boundary
    PtrList<shape> shapes_;

public:

    // Declare name of the class and its debug switch
    TypeName("immersedBoundaryMethod");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        immersedBoundaryMethod,
        dictionary,
        (
            dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from solverDict and fvMesh
    immersedBoundaryMethod
    (
        dictionary&,
        const fvMesh&,
        bool = false
    );

    // Selector from dictionary
    static autoPtr<immersedBoundaryMethod<Type,MeshType>> New
    (
        dictionary&,
        const fvMesh& fvMsh
    );

    // Destructor
    virtual ~immersedBoundaryMethod();

    // Returns true if the position vector (x,y,z)
    // is inside the immersed boundary
    bool isInside(vector);

    // Returns dimensionless wall distance given two points
    // The first point is a fluid cell and the second point is an IB cell
    scalar wallDistance(vector, vector);

    // Returns the distance from a ghost cell to the immersed boundary in
    // the wall-normal direction
    scalar wallNormalDistance(vector);

    // Returns the coordinates of the mirror point in the wall-normal direction
    // of a point located inside of the immersed boundary
    vector mirrorPoint(vector);

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
    tmp<meshField<scalar,colocated>> IBMSource
    (
        const meshField<scalar,staggered>&
    );

    // Inline access functions
    inline bool massSourceActive() {return massSourceActive_; }
    inline bool Jac() {return JacobiGhostMethod_; }
    inline meshField<label,MeshType>& mask() { return mask_; }
    inline meshField<label,MeshType>& ghostMask() { return ghostMask_; }

};

}

}

}

#ifndef solver_H
#ifdef NoRepository
    #include "immersedBoundaryMethod.C"
#endif
#endif

#endif
