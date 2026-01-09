#ifndef mappedBoundaryCondition_H
#define mappedBoundaryCondition_H

#include "mappedBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class mappedBoundaryCondition;


// Colocated

template<class Type>
class mappedBoundaryCondition<Type,colocated>
:
    public mappedBoundaryConditionBase<Type,colocated>
{
public:

    // Runtime type information
    TypeName("mapped");

    // Use constructors from base
    using mappedBoundaryConditionBase<Type,colocated>::
        mappedBoundaryConditionBase;

    // Construct from components
    mappedBoundaryCondition
    (
        const meshLevel<Type,colocated>&,
        const boundary&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new mappedBoundaryCondition<Type,colocated>(*this)
        );
    }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();
};


// Staggered

template<class Type>
class mappedBoundaryCondition<Type,staggered>
:
    public mappedBoundaryConditionBase<Type,staggered>
{
public:

    // Runtime type information
    TypeName("mapped");

    // Use constructors from base
    using mappedBoundaryConditionBase<Type,staggered>::
        mappedBoundaryConditionBase;

    // Construct from components
    mappedBoundaryCondition
    (
        const meshLevel<Type,staggered>&,
        const boundary&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new mappedBoundaryCondition<Type,staggered>(*this)
        );
    }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();
};

}

}

}

#ifdef NoRepository
    #include "mappedBoundaryCondition.C"
#endif

#endif
