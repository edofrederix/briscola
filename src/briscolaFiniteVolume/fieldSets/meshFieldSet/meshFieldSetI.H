#include "meshFieldSetReuseFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Max/min

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::max(const Type& v)
{
    forAll(*this, i)
        listType::operator[](i).max(v);
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::min(const Type& v)
{
    forAll(*this, i)
        listType::operator[](i).min(v);
}

// Member operators

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    forAll(*this, i)
        listType::operator[](i) = set[i];
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    if (tSet.isTmp() && tSet->unique())
    {
        meshFieldSet<Type,MeshType,N>& set =
            const_cast<meshFieldSet<Type,MeshType,N>&>(tSet());

        transfer(set);

        tSet.clear();
    }
    else
    {
        *this = tSet();
    }
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator=
(
    const zero
)
{
    forAll(*this, i)
        listType::operator[](i) = Zero;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator=
(
    const Type& v
)
{
    forAll(*this, i)
        listType::operator[](i) = v;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator+=
(
    const Type& v
)
{
    forAll(*this, i)
        listType::operator[](i) += v;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator+=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    forAll(*this, i)
        listType::operator[](i) += set[i];
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator+=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    *this += tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator-=
(
    const Type& v
)
{
    forAll(*this, i)
        listType::operator[](i) -= v;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator-=
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    forAll(*this, i)
        listType::operator[](i) -= set[i];
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator-=
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    *this -= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator*=
(
    const scalar s
)
{
    forAll(*this, i)
        listType::operator[](i) *= s;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator*=
(
    const meshFieldSet<scalar,MeshType,N>& set
)
{
    forAll(*this, i)
        listType::operator[](i) *= set[i];
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator*=
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet
)
{
    *this *= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator/=
(
    const scalar s
)
{
    forAll(*this, i)
        listType::operator[](i) /= s;
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator/=
(
    const meshFieldSet<scalar,MeshType,N>& set
)
{
    forAll(*this, i)
        listType::operator[](i) /= set[i];
}

template<class Type, class MeshType, direction N>
inline void meshFieldSet<Type,MeshType,N>::operator/=
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet
)
{
    *this /= tSet();

    if (tSet.isTmp())
        tSet.clear();
}

// Negation

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "-" + tSet->name()
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= -1;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

// Multiplication

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const scalar s,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + Foam::name(s) + "*" + tSet->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const meshFieldSet<scalar,MeshType,N>& s,
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const meshFieldSet<scalar,MeshType,N>& s,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + ts->name() + "*" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= ts;

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + ts->name() + "*" + tSet->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= ts;

    if (ts.isTmp())
        ts.clear();

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const meshFieldSet<Type,MeshType,N>& set,
    const meshFieldSet<scalar,MeshType,N>& s
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const meshFieldSet<scalar,MeshType,N>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "*" + s.name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const meshFieldSet<Type,MeshType,N>& set,
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "*" + ts->name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= ts;

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator*
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "*" + ts->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= ts;

    if (tSet.isTmp())
        tSet.clear();

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Scalar specialization for multiplication

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator*
(
    const meshFieldSet<scalar,MeshType,N>& set1,
    const meshFieldSet<scalar,MeshType,N>& set2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        meshFieldSet<scalar,MeshType,N>::New
        (
            "(" + set1.name() + "*" + set2.name() + ")",
            set1.fvMsh()
        );

    tRet.ref() = set1;
    tRet.ref() *= set2;

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator*
(
    const meshFieldSet<scalar,MeshType,N>& set1,
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<scalar,MeshType,N>>::New
        (
            tSet2,
            "(" + set1.name() + "*" + tSet2->name() + ")"
        );

    if (&tRet() != &tSet2())
    {
        tRet.ref() = tSet2;
    }

    tRet.ref() *= set1;

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator*
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet1,
    const meshFieldSet<scalar,MeshType,N>& set2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<scalar,MeshType,N>>::New
        (
            tSet1,
            "(" + tSet1->name() + "*" + set2.name() + ")"
        );

    if (&tRet() != &tSet1())
        tRet.ref() = tSet1;

    tRet.ref() *= set2;

    if (tSet1.isTmp())
        tSet1.clear();

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator*
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet1,
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        reuseSetTmpTmp<meshFieldSet<scalar,MeshType,N>>::New
        (
            tSet1,
            tSet2,
            "(" + tSet1->name() + "*" + tSet2->name() + ")"
        );

    if (&tRet() == &tSet1())
    {
        tRet.ref() *= tSet2;
    }
    else if (&tRet() == &tSet2())
    {
        tRet.ref() *= tSet1;
    }
    else
    {
        tRet.ref() = tSet1()*tSet2();
    }

    if (tSet1.isTmp())
        tSet1.clear();

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

// Division

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const scalar s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "/" + Foam::name(s) + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() /= s;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const scalar s,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + Foam::name(s) + "/" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = s;
    tRet.ref() /= set;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const scalar s,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + Foam::name(s) + "/" + tSet->name() + ")",
            tSet->fvMsh()
        );

    tRet.ref() = s;
    tRet.ref() /= tSet;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const meshFieldSet<Type,MeshType,N>& set,
    const meshFieldSet<scalar,MeshType,N>& s
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const meshFieldSet<scalar,MeshType,N>& s
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const meshFieldSet<Type,MeshType,N>& set,
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet(set/ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator/
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const tmp<meshFieldSet<scalar,MeshType,N>>& ts
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "/" + ts->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() /= ts;

    if (tSet.isTmp())
        tSet.clear();

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Scalar specialization for division

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator/
(
    const meshFieldSet<scalar,MeshType,N>& set1,
    const meshFieldSet<scalar,MeshType,N>& set2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        meshFieldSet<scalar,MeshType,N>::New
        (
            "(" + set1.name() + "/" + set2.name() + ")",
            set1.fvMsh()
        );

    tRet.ref() = set1;
    tRet.ref() /= set2;

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator/
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet1,
    const meshFieldSet<scalar,MeshType,N>& set2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<scalar,MeshType,N>>::New
        (
            tSet1,
            "(" + tSet1->name() + "/" + set2.name() + ")"
        );

    if (&tRet() != &tSet1())
        tRet.ref() = tSet1;

    tRet.ref() /= set2;

    if (tSet1.isTmp())
        tSet1.clear();

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator/
(
    const meshFieldSet<scalar,MeshType,N>& set1,
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet(set1/tSet2());

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

template<class MeshType, direction N>
inline tmp<meshFieldSet<scalar,MeshType,N>> operator/
(
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet1,
    const tmp<meshFieldSet<scalar,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<scalar,MeshType,N>> tRet(tSet1/tSet2());

    if (tSet1.isTmp())
        tSet1.clear();

    if (tSet2.isTmp())
        tSet2.clear();

    return tRet;
}

// Summation

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const meshFieldSet<Type,MeshType,N>& set,
    const Type& v
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "+" + Foam::name(v) + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() += v;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const Type& v
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "+" + Foam::name(v) + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() += v;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const Type& v,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + Foam::name(v) + "+" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() += v;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator+
(
    const Type& v,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + Foam::name(v) + "+" + tSet->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() += v;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

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
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
        reuseSetTmpTmp<meshFieldSet<Type,MeshType,N>>::New
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
    const meshFieldSet<Type,MeshType,N>& set,
    const Type& v
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + set.name() + "-" + Foam::name(v) + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() -= v;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet,
    const Type& v
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + tSet->name() + "-" + Foam::name(v) + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() -= v;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const Type& v,
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        meshFieldSet<Type,MeshType,N>::New
        (
            "(" + Foam::name(v) + "-" + set.name() + ")",
            set.fvMsh()
        );

    tRet.ref() = set;
    tRet.ref() *= -1;
    tRet.ref() += v;

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const Type& v,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
        (
            tSet,
            "(" + Foam::name(v) + "-" + tSet->name() + ")"
        );

    if (&tRet() != &tSet())
        tRet.ref() = tSet;

    tRet.ref() *= -1;
    tRet.ref() += v;

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

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
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const meshFieldSet<Type,MeshType,N>& set1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmp<meshFieldSet<Type,MeshType,N>>::New
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
inline tmp<meshFieldSet<Type,MeshType,N>> operator-
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet2
)
{
    tmp<meshFieldSet<Type,MeshType,N>> tRet =
        reuseSetTmpTmp<meshFieldSet<Type,MeshType,N>>::New
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

// Magnitude

#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<SCALARPRODTYPE,MeshType,N>> mag
(
    const meshFieldSet<Type,MeshType,N>& set
)
{
    tmp<meshFieldSet<SCALARPRODTYPE,MeshType,N>> tRet =
        meshFieldSet<SCALARPRODTYPE,MeshType,N>::New
        (
            "mag(" + set.name() + ")",
            set.fvMsh()
        );

    forAll(set, i)
        tRet.ref()[i] = mag(set[i]);

    return tRet;
}

template<class Type, class MeshType, direction N>
inline tmp<meshFieldSet<SCALARPRODTYPE,MeshType,N>> mag
(
    const tmp<meshFieldSet<Type,MeshType,N>>& tSet
)
{
    tmp<meshFieldSet<SCALARPRODTYPE,MeshType,N>> tRet =
        mag(tSet());

    if (tSet.isTmp())
        tSet.clear();

    return tRet;
}

#undef SCALARPRODTYPE

}

}

}
