#ifndef MittalDirichletImmersedBoundaryCondition_H
#define MittalDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"
#include "pointDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // List of coordinates of points that need to be exchanged
    // defined for each direction of the finest mesh level
    List<vectorList> exchangePoints_;

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("MittalDirichlet");

    // Construct
    MittalDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~MittalDirichletImmersedBoundaryCondition();

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
    #include "MittalDirichletImmersedBoundaryCondition.C"
#endif

#endif
