#ifndef Vreman_H
#define Vreman_H

#include "immersedBoundaryMethod.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class Vreman
:
    public immersedBoundaryMethod<Type,MeshType>
{
    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<label,MeshType> mask_;

    // The wall distance field stores the wall distances to the
    // IB in the components of the stencils
    meshField<faceScalar,MeshType> wallDist_;

    // The neighbor distance field stores the distances to neighboring
    // fluid cell for the quadratic fit on non-uniform grids
    meshField<faceScalar,MeshType> neighborDist_;

public:

    // Runtime type information
    TypeName("Vreman");

    // Construct
    Vreman
    (
        dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~Vreman();

    // Correct Jacobi points
    void correctJacobiPoints(meshLevel<Type,MeshType>&);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#ifdef NoRepository
    #include "Vreman.C"
#endif

#endif
