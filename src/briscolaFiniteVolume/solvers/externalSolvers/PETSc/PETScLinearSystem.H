#ifndef PETScLinearSystem_H
#define PETScLinearSystem_H

#include "meshDirection.H"
#include "linearSystemAggregation.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class PETScLinearSystem
:
    public linearSystemAggregation<SType,Type,MeshType>
{
    // Reference to the linear system
    const linearSystem<SType,Type,MeshType>& sys_;

    // PETSc matrices
    PtrList<Mat> matrices_;

    // Level number
    const label l_;

public:

    // Construct from components
    PETScLinearSystem
    (
        const linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label nParts
    );

    // Copy constructor
    PETScLinearSystem(const PETScLinearSystem<SType,Type,MeshType>&);

    // Destrictor
    ~PETScLinearSystem();

    // Prepare the matrix
    void prepare(const label d);

    // Return reference to the matrices

    inline Mat& matrix(const label d)
    {
        return matrices_[d];
    }

    inline const Mat& matrix(const label d) const
    {
        return matrices_[d];
    }
};

}

}

}

#ifdef NoRepository
    #include "PETScLinearSystem.C"
#endif

#endif
