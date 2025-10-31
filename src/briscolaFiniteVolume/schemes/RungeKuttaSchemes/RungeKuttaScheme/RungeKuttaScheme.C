#include "RungeKuttaScheme.H"
#include "runTimeSelectionTables.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(RungeKuttaScheme, 0);
defineRunTimeSelectionTable(RungeKuttaScheme, dictionary);

RungeKuttaScheme::RungeKuttaScheme
(
    const fvMesh& fvMsh,
    const label nStages,
    const label nSteps
)
:
    regIOobject
    (
        IOobject
        (
            "rkScheme",
            fvMsh.time().name(),
            fvMsh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    fvMsh_(fvMsh),
    dict_(fvMsh.schemeDict().subDict("rkScheme")),
    stage_(0),
    nStages_(nStages),
    nSteps_(nSteps),
    a_(nStages, nStages*nSteps, Zero),
    b_(nStages, nStages*nSteps, Zero)
{}

RungeKuttaScheme::~RungeKuttaScheme()
{}

autoPtr<RungeKuttaScheme> RungeKuttaScheme::New(const fvMesh& fvMsh)
{
    const word schemeName
    (
        fvMsh.schemeDict().subDict("rkScheme").lookup("scheme")
    );

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(schemeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Runge-Kutta scheme type "
            << schemeName << nl << nl
            << "Valid Runge-Kutta scheme types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RungeKuttaScheme>(cstrIter()(fvMsh));
}

template<class Type, class MeshType>
FastPtrList<meshField<Type,MeshType>>& RungeKuttaScheme::newStageList
(
    const word name
)
{
    typedef meshField<Type,MeshType> FieldType;

    FastPtrList<FastPtrList<FieldType>>& sources = stageSources<FieldType>();

    sources.append(new FastPtrList<FieldType>(nStages_*nSteps_));

    FastPtrList<FieldType>& list = sources.last();

    label i = 0;

    for (int step = 0; step < nSteps_; step++)
    {
        for (int stage = 0; stage < nStages_; stage++)
        {
            list.set
            (
                i,
                new FieldType
                (
                    word
                    (
                        step > 0
                      ? list[i-nStages_].name() + "_0"
                      : IOobject::groupName
                        (
                            name,
                            Foam::name(stage)
                        )
                    ),
                    fvMsh_,
                    IOobject::READ_IF_PRESENT,
                    (step < nSteps_ - 1)
                  ? IOobject::AUTO_WRITE
                  : IOobject::NO_WRITE,
                    true
                )
            );

            list[stage] = Zero;

            i++;
        }
    }

    return list;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSumA
(
    const FastPtrList<meshField<Type,MeshType>>& list
) const
{
    tmp<meshField<Type,MeshType>> tF =
        meshField<Type,MeshType>::New("K", fvMsh_);

    meshField<Type,MeshType>& F = tF.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    F = Zero;
    #endif

    for (int j = 0; j < a().n(); j++)
        if (a()(stage_-1,j) != 0.0 && (j < stage_-1 || j >= nStages_))
            F += a()(stage_-1,j)*list[j];

    return tF;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSumB
(
    const FastPtrList<meshField<Type,MeshType>>& list
) const
{
    tmp<meshField<Type,MeshType>> tF =
        meshField<Type,MeshType>::New("K", fvMsh_);

    meshField<Type,MeshType>& F = tF.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    F = Zero;
    #endif

    for (int j = 0; j < b().n(); j++)
        if (b()(stage_-1,j) != 0.0 && (j < stage_-1 || j >= nStages_))
            F += b()(stage_-1,j)*list[j];

    return tF;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSum
(
    const FastPtrList<meshField<Type,MeshType>>& listA,
    const FastPtrList<meshField<Type,MeshType>>& listB
) const
{
    return stageSumA(listA) + stageSumB(listB);
}

// Instantiate

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
                                                                               \
template FastPtrList<meshField<TYPE,MESHTYPE>>&                                \
RungeKuttaScheme::newStageList(const word);                                    \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSumA                                                    \
(const FastPtrList<meshField<TYPE,MESHTYPE>>&) const;                          \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSumB                                                    \
(const FastPtrList<meshField<TYPE,MESHTYPE>>&) const;                          \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSum                                                     \
(                                                                              \
    const FastPtrList<meshField<TYPE,MESHTYPE>>&,                              \
    const FastPtrList<meshField<TYPE,MESHTYPE>>&                               \
) const;

INSTANTIATE(scalar,colocated)
INSTANTIATE(vector,colocated)
INSTANTIATE(scalar,staggered)
INSTANTIATE(vector,staggered)

#undef INSTANTIATE

}

}

}
