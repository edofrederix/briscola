#ifndef PoissonSolvers_H
#define PoissonSolvers_H

#include "PoissonSolver.H"
#include "addToRunTimeSelectionTable.H"

#define makePoissonSolverBase(STYPE,TYPE,MESHTYPE)                             \
                                                                               \
    typedef PoissonSolver<STYPE,TYPE,MESHTYPE>                                 \
        PoissonSolver##STYPE##TYPE##MESHTYPE;                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        PoissonSolver##STYPE##TYPE##MESHTYPE,                                  \
        "PoissonSolver",                                                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        PoissonSolver##STYPE##TYPE##MESHTYPE,                                  \
        dictionary                                                             \
    );


#define makePoissonSolverType(SOLVER,STYPE,TYPE,MESHTYPE)                      \
                                                                               \
    typedef SOLVER##PoissonSolver<STYPE,TYPE,MESHTYPE>                         \
        SOLVER##PoissonSolver##STYPE##TYPE##MESHTYPE;                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SOLVER##PoissonSolver##STYPE##TYPE##MESHTYPE,                          \
        #SOLVER,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    PoissonSolver<STYPE,TYPE,MESHTYPE>::                                       \
    adddictionaryConstructorToTable                                            \
    <                                                                          \
        SOLVER##PoissonSolver##STYPE##TYPE##MESHTYPE                           \
    >                                                                          \
    add##SOLVER##PoissonSolver##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
