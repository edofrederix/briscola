#ifndef twoPhaseModels_H
#define twoPhaseModels_H

#include "twoPhaseModel.H"

#define makeTwoPhaseModel(SHORTNAME,SOLVER,VMIX,VISCOSITY,DENSITY,MESHTYPE)    \
                                                                               \
    typedef SOLVER<VMIX<VISCOSITY<DENSITY<MESHTYPE>>>>                         \
        SOLVER##VMIX##VISCOSITY##DENSITY##MESHTYPE;                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SOLVER##VMIX##VISCOSITY##DENSITY##MESHTYPE,                            \
        #SHORTNAME,                                                            \
        0                                                                      \
    );                                                                         \
                                                                               \
    twoPhaseModel::add##MESHTYPE##ConstructorToTable                           \
    <                                                                          \
        SOLVER##VMIX##VISCOSITY##DENSITY##MESHTYPE                             \
    >                                                                          \
    SOLVER##VMIX##VISCOSITY##DENSITY##MESHTYPE##ConstructorToTable_;

#endif
