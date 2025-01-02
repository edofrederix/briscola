#ifndef EigenLinearSystem_H
#define EigenLinearSystem_H

#include "EigenSolver.H"
#include "meshDirection.H"
#include "linearSystemAggregation.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class EigenLinearSystem
:
    public PtrList<EigenSolver::EigenMatrixType>,
    public linearSystemAggregation<SType,Type,MeshType>
{
    // Level number
    const label l_;

    // Diagonal bools
    List<bool> diagonal_;

public:

    // Construct from components
    EigenLinearSystem
    (
        const linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label nParts
    );

    // Copy constructor
    EigenLinearSystem(const EigenLinearSystem<SType,Type,MeshType>&);

    // Destrictor
    ~EigenLinearSystem();

    // Update the matrix
    void update(const linearSystem<SType,Type,MeshType>& sys);
};

}

}

}

#ifdef NoRepository
    #include "EigenLinearSystem.C"
#endif

#endif
