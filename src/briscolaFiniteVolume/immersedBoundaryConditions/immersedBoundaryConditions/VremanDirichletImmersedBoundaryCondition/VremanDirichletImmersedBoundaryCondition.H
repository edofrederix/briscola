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
    // List of point data exchanges defined for each level and direction
    PtrList<PtrList<cellDataExchange<MeshType>>> exchanges_;

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("VremanDirichlet");

    // Construct
    VremanDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~VremanDirichletImmersedBoundaryCondition();

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "VremanDirichletImmersedBoundaryCondition.C"
#endif

#endif
