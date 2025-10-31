#ifndef solvers_H
#define solvers_H

#include "solver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void solve(linearSystem<SType,Type,MeshType>& sys)
{
    solver<SType,Type,MeshType>::New
    (
        sys.x().name(),
        sys.x().fvMsh()
    )->solve(sys);
}

}

}

}

#define makeSolverBase(STYPE,TYPE,MESHTYPE)                                    \
                                                                               \
    typedef solver<STYPE,TYPE,MESHTYPE> solver##STYPE##TYPE##MESHTYPE;         \
                                                                               \
    typedef solver<STYPE,TYPE,MESHTYPE>::smoother                              \
    smoother##STYPE##TYPE##MESHTYPE;                                           \
                                                                               \
    typedef solver<STYPE,TYPE,MESHTYPE>::externalSolver                        \
    externalSolver##STYPE##TYPE##MESHTYPE;                                     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        solver##STYPE##TYPE##MESHTYPE,                                         \
        "solver",                                                              \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        solver##STYPE##TYPE##MESHTYPE,                                         \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        smoother##STYPE##TYPE##MESHTYPE,                                       \
        "smoother",                                                            \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        smoother##STYPE##TYPE##MESHTYPE,                                       \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        externalSolver##STYPE##TYPE##MESHTYPE,                                 \
        "smoother",                                                            \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        externalSolver##STYPE##TYPE##MESHTYPE,                                 \
        dictionary                                                             \
    );

#define makeSolverType(SOLVER,STYPE,TYPE,MESHTYPE)                             \
                                                                               \
    typedef SOLVER<STYPE,TYPE,MESHTYPE>                                        \
        SOLVER##STYPE##TYPE##MESHTYPE;                                         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SOLVER##STYPE##TYPE##MESHTYPE,                                         \
        #SOLVER,                                                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    solver<STYPE,TYPE,MESHTYPE>::                                              \
        adddictionaryConstructorToTable<SOLVER##STYPE##TYPE##MESHTYPE>         \
            add##SOLVER##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#define makeSmootherType(SMOOTHER,STYPE,TYPE,MESHTYPE)                         \
                                                                               \
    typedef SMOOTHER##Smoother<STYPE,TYPE,MESHTYPE>                            \
        SMOOTHER##Smoother##STYPE##TYPE##MESHTYPE;                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SMOOTHER##Smoother##STYPE##TYPE##MESHTYPE,                             \
        #SMOOTHER,                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    solver<STYPE,TYPE,MESHTYPE>::smoother::                                    \
    adddictionaryConstructorToTable<SMOOTHER##Smoother##STYPE##TYPE##MESHTYPE> \
        add##SMOOTHER##Smoother##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#define makeExternalSolverType(DIRECTSOLVER,STYPE,TYPE,MESHTYPE)               \
                                                                               \
    typedef DIRECTSOLVER<STYPE,TYPE,MESHTYPE>                                  \
        DIRECTSOLVER##STYPE##TYPE##MESHTYPE;                                   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        DIRECTSOLVER##STYPE##TYPE##MESHTYPE,                                   \
        #DIRECTSOLVER,                                                         \
        0                                                                      \
    );                                                                         \
                                                                               \
    solver<STYPE,TYPE,MESHTYPE>::externalSolver::                              \
        adddictionaryConstructorToTable<DIRECTSOLVER##STYPE##TYPE##MESHTYPE>   \
            add##DIRECTSOLVER##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
