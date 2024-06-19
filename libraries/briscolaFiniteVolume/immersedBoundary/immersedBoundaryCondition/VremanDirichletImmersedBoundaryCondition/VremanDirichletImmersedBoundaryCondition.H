#ifndef VremanDirichletImmersedBoundaryCondition_H
#define VremanDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"
#include "cellDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class VremanDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // List of out of bounds indices needed to correct ghost cell values
    // defined for each direction of the finest mesh level
    List<List<labelVector>> exchangePoints_;

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("VremanDirichlet");

    // Construct
    VremanDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~VremanDirichletImmersedBoundaryCondition();

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
    #include "VremanDirichletImmersedBoundaryCondition.C"
#endif

#endif
