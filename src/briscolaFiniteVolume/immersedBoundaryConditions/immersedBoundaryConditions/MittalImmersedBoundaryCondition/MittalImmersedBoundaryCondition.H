#ifndef MittalImmersedBoundaryCondition_H
#define MittalImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"
#include "pointDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{
protected:

    // List of point data exchanges defined for each level and direction
    PtrList<PtrList<pointDataExchange<MeshType>>> exchanges_;

public:

    // Constructor
    MittalImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    MittalImmersedBoundaryCondition
    (
        const MittalImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    MittalImmersedBoundaryCondition
    (
        const MittalImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~MittalImmersedBoundaryCondition();
};

}

}

}

#ifdef NoRepository
    #include "MittalImmersedBoundaryCondition.C"
#endif

#endif
