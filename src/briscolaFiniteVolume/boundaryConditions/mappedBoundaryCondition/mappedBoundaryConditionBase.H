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
public:

    // Physical type represented by Type on MeshType
    typedef typename outerProduct<Type, typename MeshType::dimType>::type
        physicalType;

protected:

    // Enable patch average enforcement
    Switch setAverage_;

    // Average to be set, defined in physical space
    physicalType average_;

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
