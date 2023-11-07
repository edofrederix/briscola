#ifndef immersedBoundaryMethodFwd_H
#define immersedBoundaryMethodFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class immersedBoundaryMethod;

typedef immersedBoundaryMethod<scalar,colocated> colocatedScalarIBM;
typedef immersedBoundaryMethod<vector,colocated> colocatedVectorIBM;
typedef immersedBoundaryMethod<scalar,staggered> staggeredScalarIBM;
typedef immersedBoundaryMethod<vector,staggered> staggeredVectorIBM;

}

}

}

#endif
