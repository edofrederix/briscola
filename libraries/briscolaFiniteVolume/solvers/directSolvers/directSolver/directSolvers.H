#ifndef directSolvers_H
#define directSolvers_H

#include "directSolver.H"

#define makeDirectSolver(STYPE,TYPE,MESHTYPE)                               \
                                                                            \
    typedef directSolver<STYPE,TYPE,MESHTYPE>                               \
        directSolver##STYPE##TYPE##MESHTYPE;                                \
                                                                            \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        directSolver##STYPE##TYPE##MESHTYPE,                                \
        "directSolver",                                                     \
        0                                                                   \
    );                                                                      \
                                                                            \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        directSolver##STYPE##TYPE##MESHTYPE,                                \
        dictionary                                                          \
    );

#define makeDirectSolverType(SOLVER,STYPE,TYPE,MESHTYPE)                    \
                                                                            \
    typedef SOLVER<STYPE,TYPE,MESHTYPE>                                     \
        SOLVER##STYPE##TYPE##MESHTYPE;                                      \
                                                                            \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        SOLVER##STYPE##TYPE##MESHTYPE,                                      \
        #SOLVER,                                                            \
        0                                                                   \
    );                                                                      \
                                                                            \
    directSolver<STYPE,TYPE,MESHTYPE>::                                     \
        adddictionaryConstructorToTable                                     \
        <SOLVER##STYPE##TYPE##MESHTYPE>                                     \
            add##SOLVER##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
