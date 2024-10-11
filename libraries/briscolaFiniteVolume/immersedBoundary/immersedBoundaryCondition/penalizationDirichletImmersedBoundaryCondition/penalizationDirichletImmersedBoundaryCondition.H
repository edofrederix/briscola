#ifndef penalizationDirichletImmersedBoundaryCondition_H
#define penalizationDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class penalizationDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("penalizationDirichlet");

    // Construct
    penalizationDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~penalizationDirichletImmersedBoundaryCondition();

    // Correct Jacobi points
    void correctJacobiPoints(meshLevel<Type,MeshType>&) const;
};

}

}

}

#ifdef NoRepository
    #include "penalizationDirichletImmersedBoundaryCondition.C"
#endif

#endif
