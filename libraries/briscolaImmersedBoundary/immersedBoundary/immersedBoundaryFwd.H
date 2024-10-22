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

template<class MeshType>
class immersedBoundary;

typedef immersedBoundary<colocated> colocatedImmersedBoundary;
typedef immersedBoundary<staggered> staggeredImmersedBoundary;

}

}

}

#endif

