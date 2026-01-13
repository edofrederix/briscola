#ifndef boundaryConditions_H
#define boundaryConditions_H

#include "boundaryCondition.H"
#include "addToRunTimeSelectionTable.H"

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
    );                                                                         \
                                                                               \
    boundaryCondition<TYPE,MESHTYPE>::                                         \
        adddictionaryConstructorToTable<BC##BoundaryCondition##TYPE##MESHTYPE> \
        add##BC##BoundaryCondition##TYPE##MESHTYPE##ConstructorToTable_;

#endif
