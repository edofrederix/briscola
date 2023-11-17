#ifndef Vreman_H
#define Vreman_H

#include "immersedBoundaryMethod.H"
#include "cellDataExchange.H"

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
    // The wall distance field stores the wall distances to the
    // IB in the components of the stencils
    meshField<faceScalar,MeshType> wallDist_;

    // The neighbor distance field stores the distances to neighboring
    // fluid cell for the quadratic fit on non-uniform grids
    meshField<faceScalar,MeshType> neighborDist_;

    // List of out of bounds indices needed to correct ghost cell values
    // defined for each level and direction
    List<List<List<labelVector>>> exchangePoints_;

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

}

}

}

#ifdef NoRepository
    #include "Vreman.C"
#endif

#endif
