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

RungeKuttaScheme::RungeKuttaScheme(const fvMesh& fvMsh)
:
    regIOobject
    (
        IOobject
        (
            "rkScheme",
            fvMsh.time().timeName(),
            fvMsh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    fvMsh_(fvMsh),
    dict_(fvMsh.schemeDict().subDict("rkScheme")),
    stage_(0)
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
PtrList<meshField<Type,MeshType>> RungeKuttaScheme::stageList
(
    const word name
) const
{
    PtrList<meshField<Type,MeshType>> list(numberOfStages());

    forAll(list, stage)
    {
        list.set
        (
            stage,
            new meshField<Type,MeshType>
            (
                IOobject::groupName
                (
                    name,
                    Foam::name(int(stage+1))
                ),
                fvMsh_
            )
        );

        list[stage] = Zero;
    }

    return list;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSumA
(
    const PtrList<meshField<Type,MeshType>>& list
) const
{
    tmp<meshField<Type,MeshType>> tF
    (
        new meshField<Type,MeshType>
        (
            a()[stage_-1][0]*list[0]
        )
    );

    meshField<Type,MeshType>& F = tF.ref();

    for (int i = 2; i < stage_; i++)
        if (a()[stage_-1][i-1] != 0.0)
            F += a()[stage_-1][i-1]*list[i-1];

    return tF;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSumB
(
    const PtrList<meshField<Type,MeshType>>& list
) const
{
    tmp<meshField<Type,MeshType>> tF
    (
        new meshField<Type,MeshType>
        (
            b()[stage_-1][0]*list[0]
        )
    );

    meshField<Type,MeshType>& F = tF.ref();

    for (int i = 2; i < stage_; i++)
        if (b()[stage_-1][i-1] != 0.0)
            F += b()[stage_-1][i-1]*list[i-1];

    return tF;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> RungeKuttaScheme::stageSum
(
    const PtrList<meshField<Type,MeshType>>& listA,
    const PtrList<meshField<Type,MeshType>>& listB
) const
{
    return stageSumA(listA) + stageSumB(listB);
}

// Instantiate

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
                                                                               \
template PtrList<meshField<TYPE,MESHTYPE>>                                     \
RungeKuttaScheme::stageList(const word) const;                                 \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSumA(const PtrList<meshField<TYPE,MESHTYPE>>&) const;   \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSumB(const PtrList<meshField<TYPE,MESHTYPE>>&) const;   \
                                                                               \
template tmp<meshField<TYPE,MESHTYPE>>                                         \
RungeKuttaScheme::stageSum                                                     \
(                                                                              \
    const PtrList<meshField<TYPE,MESHTYPE>>&,                                  \
    const PtrList<meshField<TYPE,MESHTYPE>>&                                   \
) const;

INSTANTIATE(scalar,colocated)
INSTANTIATE(vector,colocated)
INSTANTIATE(tensor,colocated)
INSTANTIATE(sphericalTensor,colocated)
INSTANTIATE(symmTensor,colocated)
INSTANTIATE(diagTensor,colocated)
INSTANTIATE(faceScalar,colocated)
INSTANTIATE(faceVector,colocated)

INSTANTIATE(scalar,staggered)
INSTANTIATE(vector,staggered)
INSTANTIATE(tensor,staggered)
INSTANTIATE(sphericalTensor,staggered)
INSTANTIATE(symmTensor,staggered)
INSTANTIATE(diagTensor,staggered)
INSTANTIATE(faceScalar,staggered)
INSTANTIATE(faceVector,staggered)

#undef INSTANTIATE

}

}

}
