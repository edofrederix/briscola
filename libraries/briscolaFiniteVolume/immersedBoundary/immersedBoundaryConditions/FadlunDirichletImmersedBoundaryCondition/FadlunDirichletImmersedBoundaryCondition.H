#ifndef FadlunDirichletImmersedBoundaryCondition_H
#define FadlunDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class FadlunDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("FadlunDirichlet");

    // Construct
    FadlunDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
    );

    // Destructor
    virtual ~FadlunDirichletImmersedBoundaryCondition();

    // Correct Jacobi points
    void correctJacobiPoints(meshLevel<Type,MeshType>&) const;
};

}

}

}

#ifdef NoRepository
    #include "FadlunDirichletImmersedBoundaryCondition.C"
#endif

#endif
