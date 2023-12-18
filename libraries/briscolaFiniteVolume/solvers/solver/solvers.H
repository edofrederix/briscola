#ifndef solvers_H
#define solvers_H

#include "solver.H"

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

#define makeSolver(STYPE,TYPE,MESHTYPE)                                        \
                                                                               \
    typedef solver<STYPE,TYPE,MESHTYPE> solver##STYPE##TYPE##MESHTYPE;         \
                                                                               \
    typedef solver<STYPE,TYPE,MESHTYPE>::smoother                              \
    smoother##STYPE##TYPE##MESHTYPE;                                           \
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
    typedef SMOOTHER<STYPE,TYPE,MESHTYPE>                                      \
        SMOOTHER##STYPE##TYPE##MESHTYPE;                                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SMOOTHER##STYPE##TYPE##MESHTYPE,                                       \
        #SMOOTHER,                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    solver<STYPE,TYPE,MESHTYPE>::smoother::                                    \
        adddictionaryConstructorToTable<SMOOTHER##STYPE##TYPE##MESHTYPE>       \
            add##SMOOTHER##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
