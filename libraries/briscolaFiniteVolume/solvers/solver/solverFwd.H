#ifndef solverFwd_H
#define solverFwd_H

#include "NamedEnum.H"
#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class solver;

typedef solver<diagStencil,scalar,colocated> colocatedScalarDiagSolver;
typedef solver<diagStencil,scalar,staggered> staggeredScalarDiagSolver;

typedef solver<diagStencil,vector,colocated> colocatedVectorDiagSolver;
typedef solver<diagStencil,vector,staggered> staggeredVectorDiagSolver;


typedef solver<stencil,scalar,colocated> colocatedScalarSolver;
typedef solver<stencil,scalar,staggered> staggeredScalarSolver;

typedef solver<stencil,vector,colocated> colocatedVectorSolver;
typedef solver<stencil,vector,staggered> staggeredVectorSolver;

}

}

}

#endif
