#ifndef boundaryConditions_H
#define boundaryConditions_H

#include "boundaryCondition.H"

#define makeBoundaryConditionBase(Type,MeshType)                                \
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

#define makeBoundaryConditionBases(MESHTYPE)                                    \
                                                                                \
    makeBoundaryConditionBase(label,MESHTYPE)                                   \
    makeBoundaryConditionBase(scalar,MESHTYPE)                                  \
    makeBoundaryConditionBase(faceScalar,MESHTYPE)                              \
    makeBoundaryConditionBase(edgeScalar,MESHTYPE)                              \
    makeBoundaryConditionBase(vertexScalar,MESHTYPE)                            \
    makeBoundaryConditionBase(vector,MESHTYPE)                                  \
    makeBoundaryConditionBase(faceVector,MESHTYPE)                              \
    makeBoundaryConditionBase(edgeVector,MESHTYPE)                              \
    makeBoundaryConditionBase(vertexVector,MESHTYPE)                            \
    makeBoundaryConditionBase(tensor,MESHTYPE)                                  \
    makeBoundaryConditionBase(sphericalTensor,MESHTYPE)                         \
    makeBoundaryConditionBase(symmTensor,MESHTYPE)                              \
    makeBoundaryConditionBase(diagTensor,MESHTYPE)                              \
    makeBoundaryConditionBase(stencil,MESHTYPE)                                 \
    makeBoundaryConditionBase(diagStencil,MESHTYPE)

#define makeBoundaryConditionTypeWithName(BC, NAME, Type, MeshType)             \
                                                                                \
    typedef BC##BoundaryCondition<Type,MeshType>                                \
        BC##BoundaryCondition##Type##MeshType;                                  \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        BC##BoundaryCondition##Type##MeshType,                                  \
        #NAME,                                                                  \
        0                                                                       \
    );                                                                          \
                                                                                \
    boundaryCondition<Type,MeshType>::                                          \
        adddictionaryConstructorToTable<BC##BoundaryCondition##Type##MeshType>  \
        add##BC##BoundaryCondition##Type##MeshType##ConstructorToTable_;

#define makeBoundaryConditionType(BC, Type, MeshType)                           \
    makeBoundaryConditionTypeWithName(BC, BC, Type, MeshType)

#define makeBoundaryConditionTypes(BC,MESHTYPE)                                 \
                                                                                \
    makeBoundaryConditionType(BC,label,MESHTYPE)                                \
    makeBoundaryConditionType(BC,scalar,MESHTYPE)                               \
    makeBoundaryConditionType(BC,faceScalar,MESHTYPE)                           \
    makeBoundaryConditionType(BC,edgeScalar,MESHTYPE)                           \
    makeBoundaryConditionType(BC,vertexScalar,MESHTYPE)                         \
    makeBoundaryConditionType(BC,vector,MESHTYPE)                               \
    makeBoundaryConditionType(BC,faceVector,MESHTYPE)                           \
    makeBoundaryConditionType(BC,edgeVector,MESHTYPE)                           \
    makeBoundaryConditionType(BC,vertexVector,MESHTYPE)                         \
    makeBoundaryConditionType(BC,tensor,MESHTYPE)                               \
    makeBoundaryConditionType(BC,sphericalTensor,MESHTYPE)                      \
    makeBoundaryConditionType(BC,symmTensor,MESHTYPE)                           \
    makeBoundaryConditionType(BC,diagTensor,MESHTYPE)                           \
    makeBoundaryConditionType(BC,stencil,MESHTYPE)                              \
    makeBoundaryConditionType(BC,diagStencil,MESHTYPE)

#endif
