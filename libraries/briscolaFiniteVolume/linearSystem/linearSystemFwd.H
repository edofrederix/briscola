#ifndef linearSystemFwd_H
#define linearSystemFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class linearSystem;

typedef linearSystem<symmStencil,scalar,colocated> colocatedScalarSymmSystem;
typedef linearSystem<symmStencil,scalar,staggered> staggeredScalarSymmSystem;

typedef linearSystem<symmStencil,vector,colocated> colocatedVectorSymmSystem;
typedef linearSystem<symmStencil,vector,staggered> staggeredVectorSymmSystem;


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
