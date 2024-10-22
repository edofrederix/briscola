#include "immersedBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

#include "emptyImmersedBoundaryCondition.H"
#include "penalizationDirichletImmersedBoundaryCondition.H"
#include "FadlunDirichletImmersedBoundaryCondition.H"
#include "MittalDirichletImmersedBoundaryCondition.H"
#include "MittalNeumannImmersedBoundaryCondition.H"
#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeIBC(label,colocated);
makeIBC(scalar,colocated);
makeIBC(faceScalar,colocated);
makeIBC(lowerFaceScalar,colocated);
makeIBC(edgeScalar,colocated);
makeIBC(vertexScalar,colocated);
makeIBC(vector,colocated);
makeIBC(faceVector,colocated);
makeIBC(lowerFaceVector,colocated);
makeIBC(edgeVector,colocated);
makeIBC(vertexVector,colocated);
makeIBC(tensor,colocated);
makeIBC(sphericalTensor,colocated)
makeIBC(symmTensor,colocated);
makeIBC(diagTensor,colocated);
makeIBC(stencil,colocated);
makeIBC(symmStencil,colocated);
makeIBC(diagStencil,colocated);

makeIBC(label,staggered);
makeIBC(scalar,staggered);
makeIBC(faceScalar,staggered);
makeIBC(lowerFaceScalar,staggered);
makeIBC(edgeScalar,staggered);
makeIBC(vertexScalar,staggered);
makeIBC(vector,staggered);
makeIBC(faceVector,staggered);
makeIBC(lowerFaceVector,staggered);
makeIBC(edgeVector,staggered);
makeIBC(vertexVector,staggered);
makeIBC(tensor,staggered);
makeIBC(sphericalTensor,staggered)
makeIBC(symmTensor,staggered);
makeIBC(diagTensor,staggered);
makeIBC(stencil,staggered);
makeIBC(symmStencil,staggered);
makeIBC(diagStencil,staggered);

makeIBCType(empty,scalar,colocated);
makeIBCType(empty,scalar,staggered);
makeIBCType(empty,vector,colocated);
makeIBCType(empty,vector,staggered);

makeIBCType(penalizationDirichlet,scalar,colocated);
makeIBCType(penalizationDirichlet,scalar,staggered);
makeIBCType(penalizationDirichlet,vector,colocated);
makeIBCType(penalizationDirichlet,vector,staggered);

makeIBCType(FadlunDirichlet,scalar,colocated);
makeIBCType(FadlunDirichlet,scalar,staggered);
makeIBCType(FadlunDirichlet,vector,colocated);
makeIBCType(FadlunDirichlet,vector,staggered);

makeIBCType(MittalDirichlet,scalar,colocated);
makeIBCType(MittalDirichlet,scalar,staggered);
makeIBCType(MittalDirichlet,vector,colocated);
makeIBCType(MittalDirichlet,vector,staggered);

makeIBCType(MittalNeumann,scalar,colocated);
makeIBCType(MittalNeumann,scalar,staggered);
makeIBCType(MittalNeumann,vector,colocated);
makeIBCType(MittalNeumann,vector,staggered);

makeIBCType(VremanDirichlet,scalar,colocated);
makeIBCType(VremanDirichlet,scalar,staggered);
makeIBCType(VremanDirichlet,vector,colocated);
makeIBCType(VremanDirichlet,vector,staggered);

}

}

}
