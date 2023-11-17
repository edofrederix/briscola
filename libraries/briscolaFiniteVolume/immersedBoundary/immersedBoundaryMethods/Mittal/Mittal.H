#ifndef Mittal_H
#define Mittal_H

#include "immersedBoundaryMethod.H"
#include "pointDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class Mittal
:
    public immersedBoundaryMethod<Type,MeshType>
{

    // List of coordinates of points that need to be exchanged
    // defined for each level and direction
    List<List<vectorList>> exchangePoints_;

    // Field containing the coordinates of the ghost mirror points
    // stored in the ghost cells
    meshField<vector,MeshType> mirrorPoints_;

public:

    // Runtime type information
    TypeName("Mittal");

    // Construct
    Mittal
    (
        dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~Mittal();

    // Correct Jacobi points
    void correctJacobiPoints(meshLevel<Type,MeshType>&);

};

}

}

}

#ifdef NoRepository
    #include "Mittal.C"
#endif

#endif
