#ifndef boundaryExchange_H
#define boundaryExchange_H

#include "meshField.H"
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
    // Reference to the field
    meshField<Type,MeshType>& field_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Pointers to parallel/periodic boundaries. Periodic boundaries are derived
    // from parallel ones.
    DynamicList<const parallelBoundaryCondition<Type,MeshType>*> boundaryPtrs_;

    // Buffers
    List<Type> sendBuffers_;
    List<Type> recvBuffers_;

    // MPI data type object
    autoPtr<pMPI<Type>> mpiDataTypePtr_;

    // Initialize buffers such that they are large enough the handle an
    // all-level boundary exchange for deep fields
    void initBuffers();

    // Correct the boundaries at given level and of given type. If the level is
    // negative then all field levels (which is one for shallow fields) are
    // corrected. If baseType is -1 then all types (periodic and parallel) are
    // corrected.
    void correct(const label l, const label baseType);

public:

    // Runtime type information
    TypeName("boundaryExchange");

    // Construct from field
    boundaryExchange(meshField<Type,MeshType>&);

    // Copy constructor
    boundaryExchange(boundaryExchange&);

    // Destructor
    virtual ~boundaryExchange();

    // Access to the field
    inline meshField<Type,MeshType>& field() { return field_; }

    // Const access to the field
    inline const meshField<Type,MeshType>& field() const { return field_; }

    // Correct both parallel and periodic boundaries of the given level
    inline void correct(const label l = 0)
    {
        this->correct(l,-1);
    }

    // Correct only parallel boundaries of the given level
    inline void correctParallel(const label l = 0)
    {
        this->correct(l,PARALLELBC);
    }

    // Correct only periodic boundaries of the given level
    inline void correctPeriodic(const label l = 0)
    {
        this->correct(l,PERIODICBC);
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
