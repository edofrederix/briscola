#ifndef twoPhaseModels_H
#define twoPhaseModels_H

#include "twoPhaseModel.H"

#define makeTwoPhaseModel(SHORTNAME,SOLVER,VISCOSITY,DENSITY)                   \
                                                                                \
    typedef SOLVER<VISCOSITY<DENSITY>>                                          \
        SOLVER##VISCOSITY##DENSITY;                                             \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SOLVER##VISCOSITY##DENSITY,                                             \
        #SHORTNAME,                                                             \
        0                                                                       \
    );                                                                          \
                                                                                \
    twoPhaseModel::adddictionaryConstructorToTable                              \
    <                                                                           \
        SOLVER##VISCOSITY##DENSITY                                              \
    >                                                                           \
    SOLVER##VISCOSITY##DENSITY##ConstructorToTable_;

#endif