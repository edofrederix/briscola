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

#define makeSolver(STYPE,TYPE,MESHTYPE)                                         \
                                                                                \
    typedef solver<STYPE,TYPE,MESHTYPE> solver##STYPE##TYPE##MESHTYPE;          \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        solver##STYPE##TYPE##MESHTYPE,                                          \
        "solver",                                                               \
        0                                                                       \
    );                                                                          \
                                                                                \
    defineTemplateRunTimeSelectionTable                                         \
    (                                                                           \
        solver##STYPE##TYPE##MESHTYPE,                                          \
        dictionary                                                              \
    );


#define makeSolverType(SOLVER,STYPE,TYPE,MESHTYPE)                              \
                                                                                \
    typedef SOLVER##Solver<STYPE,TYPE,MESHTYPE>                                 \
        SOLVER##Solver##STYPE##TYPE##MESHTYPE;                                  \
                                                                                \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        SOLVER##Solver##STYPE##TYPE##MESHTYPE,                                  \
        #SOLVER,                                                                \
        0                                                                       \
    );                                                                          \
                                                                                \
    solver<STYPE,TYPE,MESHTYPE>::                                               \
        adddictionaryConstructorToTable<SOLVER##Solver##STYPE##TYPE##MESHTYPE>  \
            add##SOLVER##Solver##STYPE##TYPE##MESHTYPE##ConstructorToTable_;

#endif
