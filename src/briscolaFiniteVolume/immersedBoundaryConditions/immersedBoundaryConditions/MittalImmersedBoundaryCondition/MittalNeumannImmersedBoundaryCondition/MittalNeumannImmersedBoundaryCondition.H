#ifndef MittalNeumannImmersedBoundaryCondition_H
#define MittalNeumannImmersedBoundaryCondition_H

#include "MittalImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalNeumannImmersedBoundaryCondition
:
    public MittalImmersedBoundaryCondition<Type,MeshType>
{
    // Boundary gradients
    List<Type> boundaryGradients_;

public:

    // Runtime type information
    TypeName("MittalNeumann");

    // Construct
    MittalNeumannImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~MittalNeumannImmersedBoundaryCondition();

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "MittalNeumannImmersedBoundaryCondition.C"
#endif

#endif
