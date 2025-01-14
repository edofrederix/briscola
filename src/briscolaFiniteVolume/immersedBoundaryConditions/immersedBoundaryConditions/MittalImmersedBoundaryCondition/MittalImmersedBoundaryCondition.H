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

    // Construct
    MittalImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& mshField,
        const immersedBoundary<MeshType>& ib
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
