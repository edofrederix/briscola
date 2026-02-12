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

    // Constructor
    emptyImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    emptyImmersedBoundaryCondition
    (
        const emptyImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    emptyImmersedBoundaryCondition
    (
        const emptyImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~emptyImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new emptyImmersedBoundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone for alternative field
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new emptyImmersedBoundaryCondition<Type,MeshType>(*this, field)
        );
    }

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
