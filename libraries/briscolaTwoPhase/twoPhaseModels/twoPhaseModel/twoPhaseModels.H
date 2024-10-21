#ifndef twoPhaseModels_H
#define twoPhaseModels_H

#include "twoPhaseModel.H"

#define makeTwoPhaseModel(SHORTNAME,SOLVER,VISCOSITY,DENSITY,MESHTYPE)         \
                                                                               \
    typedef SOLVER<VISCOSITY<DENSITY<MESHTYPE>>>                               \
        SOLVER##VISCOSITY##DENSITY##MESHTYPE;                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SOLVER##VISCOSITY##DENSITY##MESHTYPE,                                  \
        #MESHTYPE #SHORTNAME,                                                  \
        0                                                                      \
    );                                                                         \
                                                                               \
    twoPhaseModel::adddictionaryConstructorToTable                             \
    <                                                                          \
        SOLVER##VISCOSITY##DENSITY##MESHTYPE                                   \
    >                                                                          \
    SOLVER##VISCOSITY##DENSITY##MESHTYPE##ConstructorToTable_;

#endif
