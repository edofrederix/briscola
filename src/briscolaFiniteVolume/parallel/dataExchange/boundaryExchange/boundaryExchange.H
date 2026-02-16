#ifndef boundaryExchange_H
#define boundaryExchange_H

#include "meshLevel.H"
#include "parallelBoundaryCondition.H"
#include "DynamicList.H"

namespace Foam
{

template<class Type>
class pMPI;

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class boundaryExchange
{
    // Reference to the level
    meshLevel<Type,MeshType>& level_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Pointers to parallel/periodic boundaries. Nullptr entries are
    // non-parallel and can be skipped.
    List<parallelBoundaryCondition<Type,MeshType>*> boundaryPtrs_;

    // Buffers
    List<Type> sendBuffers_;
    List<Type> recvBuffers_;

    // MPI data type object
    autoPtr<pMPI<Type>> mpiDataTypePtr_;

    // Initialize buffers
    void initBuffers();

    // Correct the boundaries of given type. If baseType is -1 then all types
    // (periodic and parallel) are corrected.
    void correct(const label baseType);

public:

    // Runtime type information
    TypeName("boundaryExchange");

    // Construct from level
    boundaryExchange(meshLevel<Type,MeshType>&);

    // Copy constructor
    boundaryExchange(boundaryExchange&);

    // Destructor
    virtual ~boundaryExchange();

    // Access to the level
    inline meshLevel<Type,MeshType>& level() { return level_; }

    // Const access to the level
    inline const meshLevel<Type,MeshType>& level() const { return level_; }

    // Correct both parallel and periodic boundaries
    inline void correct()
    {
        this->correct(-1);
    }

    // Correct only parallel boundaries
    inline void correctParallel()
    {
        this->correct(PARALLELBC);
    }

    // Correct only periodic boundaries
    inline void correctPeriodic()
    {
        this->correct(PERIODICBC);
    }
};

#define makeboundaryExchange(TYPE,MESHTYPE)                                    \
                                                                               \
    typedef boundaryExchange<TYPE,MESHTYPE>                                    \
        boundaryExchange##TYPE##MESHTYPE;                                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        boundaryExchange##TYPE##MESHTYPE,                                      \
        "boundaryExchange",                                                    \
        0                                                                      \
    )                                                                          \
                                                                               \
    template class boundaryExchange<TYPE,MESHTYPE>;

#define makeboundaryExchanges(MESHTYPE)                                        \
                                                                               \
    makeboundaryExchange(label,MESHTYPE)                                       \
    makeboundaryExchange(scalar,MESHTYPE)                                      \
    makeboundaryExchange(faceScalar,MESHTYPE)                                  \
    makeboundaryExchange(edgeScalar,MESHTYPE)                                  \
    makeboundaryExchange(vertexScalar,MESHTYPE)                                \
    makeboundaryExchange(vector,MESHTYPE)                                      \
    makeboundaryExchange(faceVector,MESHTYPE)                                  \
    makeboundaryExchange(edgeVector,MESHTYPE)                                  \
    makeboundaryExchange(vertexVector,MESHTYPE)                                \
    makeboundaryExchange(tensor,MESHTYPE)                                      \
    makeboundaryExchange(sphericalTensor,MESHTYPE)                             \
    makeboundaryExchange(symmTensor,MESHTYPE)                                  \
    makeboundaryExchange(diagTensor,MESHTYPE)                                  \
    makeboundaryExchange(stencil,MESHTYPE)                                     \
    makeboundaryExchange(diagStencil,MESHTYPE)

}

}

}

#endif
