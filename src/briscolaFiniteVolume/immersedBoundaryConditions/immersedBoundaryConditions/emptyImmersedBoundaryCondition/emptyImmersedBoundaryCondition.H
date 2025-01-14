#ifndef emptyImmersedBoundaryCondition_H
#define emptyImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("empty");

    // Construct
    emptyImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~emptyImmersedBoundaryCondition();

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "emptyImmersedBoundaryCondition.C"
#endif

#endif
