#ifndef MittalNeumannImmersedBoundaryCondition_H
#define MittalNeumannImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"
#include "pointDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalNeumannImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // List of coordinates of points that need to be exchanged
    // defined for each direction of the finest mesh level
    List<vectorList> exchangePoints_;

    // Boundary gradients
    List<Type> boundaryGradients_;

public:

    // Runtime type information
    TypeName("MittalNeumann");

    // Construct
    MittalNeumannImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~MittalNeumannImmersedBoundaryCondition();

    // Correct Jacobi points
    void correctJacobiPoints
    (
        meshLevel<Type,MeshType>&
    ) const;

};

}

}

}

#ifdef NoRepository
    #include "MittalNeumannImmersedBoundaryCondition.C"
#endif

#endif
