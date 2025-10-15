#ifndef meshFieldSetFwd_H
#define meshFieldSetFwd_H

#include "direction.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType, direction N>
class meshFieldSet;


// typedef meshFieldSet<diagStencil,scalar,colocated> colocatedScalarDiagSystem;
// typedef meshFieldSet<diagStencil,scalar,staggered> staggeredScalarDiagSystem;

// typedef meshFieldSet<diagStencil,vector,colocated> colocatedVectorDiagSystem;
// typedef meshFieldSet<diagStencil,vector,staggered> staggeredVectorDiagSystem;


// typedef meshFieldSet<stencil,scalar,colocated> colocatedScalarSystem;
// typedef meshFieldSet<stencil,scalar,staggered> staggeredScalarSystem;

// typedef meshFieldSet<stencil,vector,colocated> colocatedVectorSystem;
// typedef meshFieldSet<stencil,vector,staggered> staggeredVectorSystem;

}

}

}

#endif
