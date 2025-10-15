#include "meshFieldSetReuseFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Negation

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator-
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "-" + set.name(),
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= -1;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "-" + tSet->name()
        );

    for (direction s = 0; s < N; s++)
        tRet->operator[](s) = tSet()[s];

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

// Multiplication

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const scalar s,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + Foam::name(s) + "*" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const scalar s,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + Foam::name(s) + "*" + tSet->name() + ")"
        );

    for (direction s = 0; s < N; s++)
        tRet->operator[](s) = s*tSet()[s];

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const meshFieldSet<Type,MeshType,N>& set,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "*" + Foam::name(s) + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + tSet->name() + "*" + Foam::name(s) + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const meshField<scalar,MeshType>& s,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + s.name() + "*" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const meshField<scalar,MeshType>& s,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + s.name() + "*" + tSet->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const tmp<meshField<scalar,MeshType>>& ts,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(ts()*set);

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const tmp<meshField<scalar,MeshType>>& ts,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(ts()*tSet);

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const meshFieldSet<Type,MeshType,N>& set,
    const meshField<scalar,MeshType>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "*" + s.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const meshField<scalar,MeshType>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + tSet->name() + "*" + s.name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const meshFieldSet<Type,MeshType,N>& set,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(set*ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(tSet*ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Division

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const meshFieldSet<Type,MeshType,N>& set,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "/" + Foam::name(s) + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() /= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + tSet->name() + "/" + Foam::name(s) + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() /= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const meshFieldSet<Type,MeshType,N>& set,
    const meshField<scalar,MeshType>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "/" + s.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() /= s;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const meshField<scalar,MeshType>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet,
            "(" + tSet->name() + "/" + s.name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() /= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const meshFieldSet<Type,MeshType,N>& set,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(set/ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(tSet/ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Summation

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const meshFieldSet<Type,MeshType,N>& set1,
    const meshFieldSet<Type,MeshType,N>& set2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set1.name() + "+" + set2.name() + ")",
            set1.fvMsh()
        );

    tRet.ref() = set1;
    tRet.ref() += set2;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
    const meshFieldSet<Type,MeshType,N>& set2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet1,
            "(" + tSet1->name() + "+" + set2.name() + ")"
        );

    if (&tRet() != &tSet1())
        tRet.ref() = tSet1;

    tRet.ref() += set2;

    if (tSet1.isTmp())
        tSet1.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const meshFieldSet<Type,MeshType,N>& set1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet2,
            "(" + set1.name() + "+" + tSet2->name() + ")"
        );

    if (&tRet() != & tSet2())
        tRet.ref() = tSet2;

    tRet.ref() += set1;

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmpTmp<Type,MeshType,N>::New
        (
            tSet1,
            tSet2,
            "(" + tSet1->name() + "+" + tSet2->name() + ")"
        );

    if (&tRet() == &tSet1())
    {
        tRet.ref() += tSet2;
    }
    else if(&tRet() == &tSet2())
    {
        tRet.ref() += tSet1;
    }
    else
    {
        tRet.ref() = tSet1;
        tRet.ref() += tSet2;
    }

    if (tSet1.isTmp())
        tSet1.clear();

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

// Subtraction

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const meshFieldSet<Type,MeshType,N>& set1,
    const meshFieldSet<Type,MeshType,N>& set2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set1.name() + "-" + set2.name() + ")",
            set1.fvMsh()
        );

    tRet.ref() = set1;
    tRet.ref() -= set2;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
    const meshFieldSet<Type,MeshType,N>& set2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet1,
            "(" + tSet1->name() + "-" + set2.name() + ")"
        );

    if (&tRet() != & tSet1())
        tRet.ref() = tSet1;

    tRet.ref() -= set2;

    if (tSet1.isTmp())
        tSet1.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator-
(
    const meshFieldSet<Type,MeshType,N>& set1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<Type,MeshType,N>::New
        (
            tSet2,
            "(" + set1.name() + "-" + tSet2->name() + ")"
        );

    if (&tRet() != & tSet2())
        tRet.ref() = tSet2;

    tRet.ref() *= -1;
    tRet.ref() += set1;

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>>
operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmpTmp<Type,MeshType,N>::New
        (
            tSet1,
            tSet2,
            "(" + tSet1->name() + "-" + tSet2->name() + ")"
        );

    if (&tRet() == &tSet1())
    {
        tRet.ref() -= tSet2;
    }
    else if(&tRet() == &tSet2())
    {
        tRet.ref() *= -1;
        tRet.ref() += tSet1;
    }
    else
    {
        tRet.ref() = tSet1;
        tRet.ref() -= tSet2;
    }

    if (tSet1.isTmp())
        tSet1.clear();

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

}

}

}
