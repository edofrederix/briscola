#ifndef MittalDirichletImmersedBoundaryCondition_H
#define MittalDirichletImmersedBoundaryCondition_H

#include "MittalImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalDirichletImmersedBoundaryCondition
:
    public MittalImmersedBoundaryCondition<Type,MeshType>
{
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

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "MittalDirichletImmersedBoundaryCondition.C"
#endif

#endif
