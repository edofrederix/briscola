#ifndef boundaryConditions_H
#define boundaryConditions_H

#include "boundaryCondition.H"

#define makeBoundaryCondition(Type,MeshType)                                    \
                                                                                \
    typedef boundaryCondition<Type,MeshType>                                    \
        boundaryCondition##Type##MeshType;                                      \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        boundaryCondition##Type##MeshType,                                      \
        "boundaryCondition",                                                    \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        boundaryCondition##Type##MeshType,                                      \
        dictionary                                                              \
    );

#define makeBoundaryConditionType(BC, Type, MeshType)                           \
                                                                                \
    typedef BC##BoundaryCondition<Type,MeshType>                                \
        BC##BoundaryCondition##Type##MeshType;                                  \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        BC##BoundaryCondition##Type##MeshType,                                  \
        #BC,                                                                    \
        0                                                                       \
    );                                                                          \
                                                                                \
    boundaryCondition<Type,MeshType>::                                          \
        adddictionaryConstructorToTable<BC##BoundaryCondition##Type##MeshType>  \
        add##BC##BoundaryCondition##Type##MeshType##ConstructorToTable_;

#define makeBoundaryConditionTypes(BC)                                          \
                                                                                \
    makeBoundaryConditionType(BC,label,colocated)                               \
    makeBoundaryConditionType(BC,label,staggered)                               \
                                                                                \
    makeBoundaryConditionType(BC,scalar,colocated)                              \
    makeBoundaryConditionType(BC,scalar,staggered)                              \
                                                                                \
    makeBoundaryConditionType(BC,faceScalar,colocated)                          \
    makeBoundaryConditionType(BC,faceScalar,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,edgeScalar,colocated)                          \
    makeBoundaryConditionType(BC,edgeScalar,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,vertexScalar,colocated)                        \
    makeBoundaryConditionType(BC,vertexScalar,staggered)                        \
                                                                                \
    makeBoundaryConditionType(BC,vector,colocated)                              \
    makeBoundaryConditionType(BC,vector,staggered)                              \
                                                                                \
    makeBoundaryConditionType(BC,faceVector,colocated)                          \
    makeBoundaryConditionType(BC,faceVector,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,edgeVector,colocated)                          \
    makeBoundaryConditionType(BC,edgeVector,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,vertexVector,colocated)                        \
    makeBoundaryConditionType(BC,vertexVector,staggered)                        \
                                                                                \
    makeBoundaryConditionType(BC,tensor,colocated)                              \
    makeBoundaryConditionType(BC,tensor,staggered)                              \
                                                                                \
    makeBoundaryConditionType(BC,sphericalTensor,colocated)                     \
    makeBoundaryConditionType(BC,sphericalTensor,staggered)                     \
                                                                                \
    makeBoundaryConditionType(BC,symmTensor,colocated)                          \
    makeBoundaryConditionType(BC,symmTensor,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,diagTensor,colocated)                          \
    makeBoundaryConditionType(BC,diagTensor,staggered)                          \
                                                                                \
    makeBoundaryConditionType(BC,stencil,colocated)                             \
    makeBoundaryConditionType(BC,stencil,staggered)                             \
                                                                                \
    makeBoundaryConditionType(BC,diagStencil,colocated)                         \
    makeBoundaryConditionType(BC,diagStencil,staggered)

#endif
