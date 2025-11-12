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

    // Append a list that stores all but the last stage in each step

    sources.append(new FastPtrList<FieldType>((nStages_ - 1)*nSteps_));

    FastPtrList<FieldType>& list = sources.last();

    label i = 0;

    for (int step = 0; step < nSteps_; step++)
    {
        for (int stage = 0; stage < nStages_-1; stage++)
        {
            list.set
            (
                i,
                new FieldType
                (
                    word
                    (
                        step > 0
                      ? list[i-(nStages_-1)].name() + "_0"
                      : IOobject::groupName
                        (
                            "RK:" + name,
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

            list[i] = Zero;

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

    F = Zero;

    label j = 0;

    for (int step = 0; step < nSteps_; step++)
    {
        for (int stage = 0; stage < nStages_-1; stage++)
        {
            const label k = stage + step*nStages_;

            if (a()(stage_-1,k) != 0.0)
                if (k < stage_-1 || step > 0)
                    F += a()(stage_-1,k)*list[j];

            j++;
        }
    }

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

    F = Zero;

    label j = 0;

    for (int step = 0; step < nSteps_; step++)
    {
        for (int stage = 0; stage < nStages_-1; stage++)
        {
            const label k = stage + step*nStages_;

            if (b()(stage_-1,k) != 0.0)
                if (k < stage_-1 || step > 0)
                    F += b()(stage_-1,k)*list[j];

            j++;
        }
    }

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
