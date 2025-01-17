#ifndef linearSystemAggregationFwd_H
#define linearSystemAggregationFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class linearSystemAggregation;


typedef linearSystemAggregation<diagStencil,scalar,colocated>
    colocatedScalarDiagSystemAggregation;
typedef linearSystemAggregation<diagStencil,scalar,staggered>
    staggeredScalarDiagSystemAggregation;

typedef linearSystemAggregation<diagStencil,vector,colocated>
    colocatedVectorDiagSystemAggregation;
typedef linearSystemAggregation<diagStencil,vector,staggered>
    staggeredVectorDiagSystemAggregation;


typedef linearSystemAggregation<symmStencil,scalar,colocated>
    colocatedScalarSymmSystemAggregation;
typedef linearSystemAggregation<symmStencil,scalar,staggered>
    staggeredScalarSymmSystemAggregation;

typedef linearSystemAggregation<symmStencil,vector,colocated>
    colocatedVectorSymmSystemAggregation;
typedef linearSystemAggregation<symmStencil,vector,staggered>
    staggeredVectorSymmSystemAggregation;


typedef linearSystemAggregation<stencil,scalar,colocated>
    colocatedScalarSystemAggregation;
typedef linearSystemAggregation<stencil,scalar,staggered>
    staggeredScalarSystemAggregation;

typedef linearSystemAggregation<stencil,vector,colocated>
    colocatedVectorSystemAggregation;
typedef linearSystemAggregation<stencil,vector,staggered>
    staggeredVectorSystemAggregation;

}

}

}

#endif
