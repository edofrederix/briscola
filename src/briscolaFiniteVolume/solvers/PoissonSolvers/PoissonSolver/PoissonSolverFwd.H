#ifndef PoissonSolverFwd_H
#define PoissonSolverFwd_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class PoissonSolver;

typedef PoissonSolver<diagStencil,scalar,colocated>
    colocatedScalarDiagPoissonSolver;
typedef PoissonSolver<diagStencil,scalar,staggered>
    staggeredScalarDiagPoissonSolver;

typedef PoissonSolver<diagStencil,vector,colocated>
    colocatedVectorDiagPoissonSolver;
typedef PoissonSolver<diagStencil,vector,staggered>
    staggeredVectorDiagPoissonSolver;


typedef PoissonSolver<stencil,scalar,colocated>
    colocatedScalarPoissonSolver;
typedef PoissonSolver<stencil,scalar,staggered>
    staggeredScalarPoissonSolver;

typedef PoissonSolver<stencil,vector,colocated>
    colocatedVectorPoissonSolver;
typedef PoissonSolver<stencil,vector,staggered>
    staggeredVectorPoissonSolver;

}

}

}

#endif
