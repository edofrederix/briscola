#ifndef mappedBoundaryConditionBase_H
#define mappedBoundaryConditionBase_H

#include "DirichletBoundaryCondition.H"
#include "pointDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class mappedBoundaryConditionBase
:
    public DirichletBoundaryCondition<Type,MeshType>
{
protected:

    // List of switches to enable patch average enforcement
    List<Switch> setAverages_;

    // List of mean values to be set
    List<Type> averages_;

    // Point data exchanges
    PtrList<pointDataExchange<MeshType>> exchanges_;

public:

    // Use constructors from base
    using DirichletBoundaryCondition<Type,MeshType>::
        DirichletBoundaryCondition;

    // Construct from components
    mappedBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );
};

}

}

}

#ifdef NoRepository
    #include "mappedBoundaryConditionBase.C"
#endif

#endif
