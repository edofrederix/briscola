#ifndef linearSystems_H
#define linearSystems_H

#include "linearSystem.H"
#include "colocatedFields.H"
#include "staggeredFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef linearSystem<diagStencil,scalar,colocated> colocatedScalarDiagSystem;
typedef linearSystem<diagStencil,scalar,staggered> staggeredScalarDiagSystem;

typedef linearSystem<diagStencil,vector,colocated> colocatedVectorDiagSystem;
typedef linearSystem<diagStencil,vector,staggered> staggeredVectorDiagSystem;

typedef linearSystem<stencil,scalar,colocated> colocatedScalarSystem;
typedef linearSystem<stencil,scalar,staggered> staggeredScalarSystem;

typedef linearSystem<stencil,vector,colocated> colocatedVectorSystem;
typedef linearSystem<stencil,vector,staggered> staggeredVectorSystem;

}

}

}

#endif
