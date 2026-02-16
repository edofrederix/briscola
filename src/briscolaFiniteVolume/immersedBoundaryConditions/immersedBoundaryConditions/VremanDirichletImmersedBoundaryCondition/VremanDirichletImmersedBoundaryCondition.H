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

    // Constructor
    VremanDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    VremanDirichletImmersedBoundaryCondition
    (
        const VremanDirichletImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    VremanDirichletImmersedBoundaryCondition
    (
        const VremanDirichletImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~VremanDirichletImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new VremanDirichletImmersedBoundaryCondition<Type,MeshType>(*this)
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
            new VremanDirichletImmersedBoundaryCondition<Type,MeshType>
            (
                *this,
                field
            )
        );
    }

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
