#ifndef TaylorBubbleInletBoundaryCondition_H
#define TaylorBubbleInletBoundaryCondition_H

#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class MeshType>
class TaylorBubbleInletBoundaryCondition;


// Colocated

template<>
class TaylorBubbleInletBoundaryCondition<colocated>
:
    public DirichletBoundaryCondition<vector,colocated>
{
public:

    // Runtime type information
    TypeName("TaylorBubbleInlet");

    // Use base class constructors
    using DirichletBoundaryCondition<vector,colocated>::
        DirichletBoundaryCondition;

    // Construct from components
    TaylorBubbleInletBoundaryCondition
    (
        const meshLevel<vector,colocated>& level,
        const boundary& b
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<vector,colocated>> clone() const
    {
        return tmp<boundaryCondition<vector,colocated>>
        (
            new TaylorBubbleInletBoundaryCondition<colocated>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<vector,colocated>> clone
    (
        const meshLevel<vector,colocated>& level
    ) const
    {
        return tmp<boundaryCondition<vector,colocated>>
        (
            new TaylorBubbleInletBoundaryCondition<colocated>
            (
                *this,
                level
            )
        );
    }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();
};


// Staggered

template<>
class TaylorBubbleInletBoundaryCondition<staggered>
:
    public DirichletBoundaryCondition<scalar,staggered>
{
public:

    // Runtime type information
    TypeName("TaylorBubbleInlet");

    // Use base class constructors
    using DirichletBoundaryCondition<scalar,staggered>::
        DirichletBoundaryCondition;

    // Construct from components
    TaylorBubbleInletBoundaryCondition
    (
        const meshLevel<scalar,staggered>& level,
        const boundary& b
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<scalar,staggered>> clone() const
    {
        return tmp<boundaryCondition<scalar,staggered>>
        (
            new TaylorBubbleInletBoundaryCondition<staggered>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<scalar,staggered>> clone
    (
        const meshLevel<scalar,staggered>& level
    ) const
    {
        return tmp<boundaryCondition<scalar,staggered>>
        (
            new TaylorBubbleInletBoundaryCondition<staggered>
            (
                *this,
                level
            )
        );
    }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();
};

}

}

}

#endif
