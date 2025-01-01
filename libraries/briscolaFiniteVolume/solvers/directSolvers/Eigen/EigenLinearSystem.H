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
public:

    // Construct from components
    EigenLinearSystem
    (
        const meshLevel<SType,MeshType>& A,
        const label nParts,
        const List<bool>& singular
    );

    // Copy constructor
    EigenLinearSystem(const EigenLinearSystem<SType,Type,MeshType>&);

    // Destrictor
    ~EigenLinearSystem();

    // Update the matrix
    void update
    (
        const meshLevel<SType,MeshType>& A,
        const List<bool>& singular
    );
};

}

}

}

#ifdef NoRepository
    #include "EigenLinearSystem.C"
#endif

#endif
