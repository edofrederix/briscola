#ifndef boundaryConditions_H
#define boundaryConditions_H

#include "boundaryCondition.H"

#define makeBoundaryConditionBase(TYPE,MESHTYPE)                               \
                                                                               \
    typedef boundaryCondition<TYPE,MESHTYPE>                                   \
        boundaryCondition##TYPE##MESHTYPE;                                     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        boundaryCondition##TYPE##MESHTYPE,                                     \
        "boundaryCondition",                                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        boundaryCondition##TYPE##MESHTYPE,                                     \
        dictionary                                                             \
    );

#define makeBoundaryConditionBases(MESHTYPE)                                   \
                                                                               \
    makeBoundaryConditionBase(label,MESHTYPE)                                  \
    makeBoundaryConditionBase(scalar,MESHTYPE)                                 \
    makeBoundaryConditionBase(faceScalar,MESHTYPE)                             \
    makeBoundaryConditionBase(lowerFaceScalar,MESHTYPE)                        \
    makeBoundaryConditionBase(edgeScalar,MESHTYPE)                             \
    makeBoundaryConditionBase(vertexScalar,MESHTYPE)                           \
    makeBoundaryConditionBase(vector,MESHTYPE)                                 \
    makeBoundaryConditionBase(faceVector,MESHTYPE)                             \
    makeBoundaryConditionBase(lowerFaceVector,MESHTYPE)                        \
    makeBoundaryConditionBase(edgeVector,MESHTYPE)                             \
    makeBoundaryConditionBase(vertexVector,MESHTYPE)                           \
    makeBoundaryConditionBase(tensor,MESHTYPE)                                 \
    makeBoundaryConditionBase(sphericalTensor,MESHTYPE)                        \
    makeBoundaryConditionBase(symmTensor,MESHTYPE)                             \
    makeBoundaryConditionBase(diagTensor,MESHTYPE)                             \
    makeBoundaryConditionBase(stencil,MESHTYPE)                                \
    makeBoundaryConditionBase(symmStencil,MESHTYPE)                            \
    makeBoundaryConditionBase(diagStencil,MESHTYPE)

#define makeBoundaryConditionType(BC,NAME,TYPE,MESHTYPE)                       \
                                                                               \
    typedef BC##BoundaryCondition<TYPE,MESHTYPE>                               \
        BC##BoundaryCondition##TYPE##MESHTYPE;                                 \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BC##BoundaryCondition##TYPE##MESHTYPE,                                 \
        #NAME,                                                                 \
        0                                                                      \
    );

#define addBoundaryConditionType(BC,TYPE,MESHTYPE)                             \
                                                                               \
    boundaryCondition<TYPE,MESHTYPE>::                                         \
        adddictionaryConstructorToTable<BC##BoundaryCondition##TYPE##MESHTYPE> \
        add##BC##BoundaryCondition##TYPE##MESHTYPE##ConstructorToTable_;

#define makeBoundaryConditionMeshedType(BC,NAME,TYPE,MESHTYPE)                 \
                                                                               \
    typedef BC##BoundaryCondition<TYPE>                                        \
        BC##BoundaryCondition##TYPE;                                           \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BC##BoundaryCondition##TYPE,                                           \
        #NAME,                                                                 \
        0                                                                      \
    );

#define addBoundaryConditionMeshedType(BC,TYPE,MESHTYPE)                       \
                                                                               \
    boundaryCondition<TYPE,MESHTYPE>::                                         \
        adddictionaryConstructorToTable<BC##BoundaryCondition##TYPE>           \
        add##BC##BoundaryCondition##TYPE##ConstructorToTable_;

#define makeBoundaryConditionTypes(BC,NAME,MESHTYPE)                           \
                                                                               \
    makeBoundaryConditionType(BC,NAME,label,MESHTYPE)                          \
    makeBoundaryConditionType(BC,NAME,scalar,MESHTYPE)                         \
    makeBoundaryConditionType(BC,NAME,faceScalar,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,lowerFaceScalar,MESHTYPE)                \
    makeBoundaryConditionType(BC,NAME,edgeScalar,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,vertexScalar,MESHTYPE)                   \
    makeBoundaryConditionType(BC,NAME,vector,MESHTYPE)                         \
    makeBoundaryConditionType(BC,NAME,faceVector,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,lowerFaceVector,MESHTYPE)                \
    makeBoundaryConditionType(BC,NAME,edgeVector,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,vertexVector,MESHTYPE)                   \
    makeBoundaryConditionType(BC,NAME,tensor,MESHTYPE)                         \
    makeBoundaryConditionType(BC,NAME,sphericalTensor,MESHTYPE)                \
    makeBoundaryConditionType(BC,NAME,symmTensor,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,diagTensor,MESHTYPE)                     \
    makeBoundaryConditionType(BC,NAME,stencil,MESHTYPE)                        \
    makeBoundaryConditionType(BC,NAME,symmStencil,MESHTYPE)                    \
    makeBoundaryConditionType(BC,NAME,diagStencil,MESHTYPE)

#define addBoundaryConditionTypes(BC,MESHTYPE)                                 \
                                                                               \
    addBoundaryConditionType(BC,label,MESHTYPE)                                \
    addBoundaryConditionType(BC,scalar,MESHTYPE)                               \
    addBoundaryConditionType(BC,faceScalar,MESHTYPE)                           \
    addBoundaryConditionType(BC,lowerFaceScalar,MESHTYPE)                      \
    addBoundaryConditionType(BC,edgeScalar,MESHTYPE)                           \
    addBoundaryConditionType(BC,vertexScalar,MESHTYPE)                         \
    addBoundaryConditionType(BC,vector,MESHTYPE)                               \
    addBoundaryConditionType(BC,faceVector,MESHTYPE)                           \
    addBoundaryConditionType(BC,lowerFaceVector,MESHTYPE)                      \
    addBoundaryConditionType(BC,edgeVector,MESHTYPE)                           \
    addBoundaryConditionType(BC,vertexVector,MESHTYPE)                         \
    addBoundaryConditionType(BC,tensor,MESHTYPE)                               \
    addBoundaryConditionType(BC,sphericalTensor,MESHTYPE)                      \
    addBoundaryConditionType(BC,symmTensor,MESHTYPE)                           \
    addBoundaryConditionType(BC,diagTensor,MESHTYPE)                           \
    addBoundaryConditionType(BC,stencil,MESHTYPE)                              \
    addBoundaryConditionType(BC,symmStencil,MESHTYPE)                          \
    addBoundaryConditionType(BC,diagStencil,MESHTYPE)

#define makeBoundaryConditionMeshedTypes(BC,NAME,MESHTYPE)                     \
                                                                               \
    makeBoundaryConditionMeshedType(BC,NAME,label,MESHTYPE)                    \
    makeBoundaryConditionMeshedType(BC,NAME,scalar,MESHTYPE)                   \
    makeBoundaryConditionMeshedType(BC,NAME,faceScalar,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,lowerFaceScalar,MESHTYPE)          \
    makeBoundaryConditionMeshedType(BC,NAME,edgeScalar,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,vertexScalar,MESHTYPE)             \
    makeBoundaryConditionMeshedType(BC,NAME,vector,MESHTYPE)                   \
    makeBoundaryConditionMeshedType(BC,NAME,faceVector,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,lowerFaceVector,MESHTYPE)          \
    makeBoundaryConditionMeshedType(BC,NAME,edgeVector,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,vertexVector,MESHTYPE)             \
    makeBoundaryConditionMeshedType(BC,NAME,tensor,MESHTYPE)                   \
    makeBoundaryConditionMeshedType(BC,NAME,sphericalTensor,MESHTYPE)          \
    makeBoundaryConditionMeshedType(BC,NAME,symmTensor,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,diagTensor,MESHTYPE)               \
    makeBoundaryConditionMeshedType(BC,NAME,stencil,MESHTYPE)                  \
    makeBoundaryConditionMeshedType(BC,NAME,symmStencil,MESHTYPE)              \
    makeBoundaryConditionMeshedType(BC,NAME,diagStencil,MESHTYPE)

#define addBoundaryConditionMeshedTypes(BC,MESHTYPE)                           \
                                                                               \
    addBoundaryConditionMeshedType(BC,label,MESHTYPE)                          \
    addBoundaryConditionMeshedType(BC,scalar,MESHTYPE)                         \
    addBoundaryConditionMeshedType(BC,faceScalar,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,lowerFaceScalar,MESHTYPE)                \
    addBoundaryConditionMeshedType(BC,edgeScalar,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,vertexScalar,MESHTYPE)                   \
    addBoundaryConditionMeshedType(BC,vector,MESHTYPE)                         \
    addBoundaryConditionMeshedType(BC,faceVector,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,lowerFaceVector,MESHTYPE)                \
    addBoundaryConditionMeshedType(BC,edgeVector,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,vertexVector,MESHTYPE)                   \
    addBoundaryConditionMeshedType(BC,tensor,MESHTYPE)                         \
    addBoundaryConditionMeshedType(BC,sphericalTensor,MESHTYPE)                \
    addBoundaryConditionMeshedType(BC,symmTensor,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,diagTensor,MESHTYPE)                     \
    addBoundaryConditionMeshedType(BC,stencil,MESHTYPE)                        \
    addBoundaryConditionMeshedType(BC,symmStencil,MESHTYPE)                    \
    addBoundaryConditionMeshedType(BC,diagStencil,MESHTYPE)

#endif
