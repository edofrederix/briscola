#ifndef immersedBoundaryConditions_H
#define immersedBoundaryConditions_H

#include "immersedBoundaryCondition.H"

#define makeIBC(TYPE,MESHTYPE)                                                  \
                                                                                \
    typedef immersedBoundaryCondition<TYPE,MESHTYPE>                            \
        immersedBoundaryCondition##TYPE##MESHTYPE;                              \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        immersedBoundaryCondition##TYPE##MESHTYPE,                              \
        "immersedBoundaryCondition",                                            \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        immersedBoundaryCondition##TYPE##MESHTYPE,                              \
        dictionary                                                              \
    );


#define makeIBCType(IBC,TYPE,MESHTYPE)                                          \
                                                                                \
    typedef IBC##ImmersedBoundaryCondition<TYPE,MESHTYPE>                       \
        IBC##TYPE##MESHTYPE;                                                    \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        IBC##TYPE##MESHTYPE,                                                    \
        #IBC,                                                                   \
        0                                                                       \
    );                                                                          \
                                                                                \
    immersedBoundaryCondition<TYPE,MESHTYPE>::                                  \
        adddictionaryConstructorToTable<IBC##TYPE##MESHTYPE>                    \
            add##IBC##TYPE##MESHTYPE##ConstructorToTable_;

#endif
