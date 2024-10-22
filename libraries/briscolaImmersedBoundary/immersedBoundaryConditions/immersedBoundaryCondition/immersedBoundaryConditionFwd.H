#ifndef immersedBoundaryConditionFwd_H
#define immersedBoundaryConditionFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class immersedBoundaryCondition;

typedef immersedBoundaryCondition<scalar,colocated> colocatedScalarIBC;
typedef immersedBoundaryCondition<vector,colocated> colocatedVectorIBC;
typedef immersedBoundaryCondition<scalar,staggered> staggeredScalarIBC;
typedef immersedBoundaryCondition<vector,staggered> staggeredVectorIBC;

}

}

}

#endif
