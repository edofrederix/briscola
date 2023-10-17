#ifndef immersedBoundaryFwd_H
#define immersedBoundaryFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class immersedBoundary;

typedef immersedBoundary<scalar,colocated> colocatedScalarImmersedBoundary;
typedef immersedBoundary<vector,colocated> colocatedVectorImmersedBoundary;
typedef immersedBoundary<scalar,staggered> staggeredScalarImmersedBoundary;
typedef immersedBoundary<vector,staggered> staggeredVectorImmersedBoundary;

}

}

}

#endif
