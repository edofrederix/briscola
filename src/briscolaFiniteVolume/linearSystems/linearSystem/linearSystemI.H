#include "linearSystemReuseFunctions.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Negation

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "-" + sys.name(),
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() *= -1;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "-" + tSys->name()
        );

    tRet->A() = -tSys->A();
    tRet->b() = -tSys->b();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Multiplication

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const scalar s,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + Foam::name(s) + "*" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() *= s;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const scalar s,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + Foam::name(s) + "*" + tSys->name() + ")"
        );

    tRet->A() = s*tSys->A();
    tRet->b() = s*tSys->b();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const linearSystem<SType,Type,MeshType>& sys,
    const scalar s
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "*" + Foam::name(s) + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() *= s;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const scalar s
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "*" + Foam::name(s) + ")"
        );

    tRet->A() = tRet->A()*s;
    tRet->b() = tRet->b()*s;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const meshField<scalar,MeshType>& s,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + s.name() + "*" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const meshField<scalar,MeshType>& s,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + s.name() + "*" + tSys->name() + ")"
        );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation faults.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    const meshField<SType,MeshType>& As = tSys->A();
    const meshField<Type,MeshType>& bs = tSys->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) = As(l,d,i,j,k)*s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) = bs(l,d,i,j,k)*s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const tmp<meshField<scalar,MeshType>>& ts,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(ts()*sys);

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const tmp<meshField<scalar,MeshType>>& ts,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(ts()*tSys);

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const linearSystem<SType,Type,MeshType>& sys,
    const meshField<scalar,MeshType>& s
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "*" + s.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const meshField<scalar,MeshType>& s
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "*" + s.name() + ")"
        );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    const meshField<SType,MeshType>& As = tSys->A();
    const meshField<Type,MeshType>& bs = tSys->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) = As(l,d,i,j,k)*s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) = bs(l,d,i,j,k)*s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const linearSystem<SType,Type,MeshType>& sys,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(sys*ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(tSys*ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Division

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const linearSystem<SType,Type,MeshType>& sys,
    const scalar s
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "/" + Foam::name(s) + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;

    tRet->A() /= s;
    tRet->b() /= s;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const scalar s
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "/" + Foam::name(s) + ")"
        );

    tRet->A() = tSys->A()/s;
    tRet->b() = tSys->b()/s;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const linearSystem<SType,Type,MeshType>& sys,
    const meshField<scalar,MeshType>& s
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "/" + s.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) /= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) /= s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const meshField<scalar,MeshType>& s
)
{
    const bool shallow = s.shallow();

    if (shallow)
        restrict(s);

    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "/" + s.name() + ")"
        );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    const meshField<SType,MeshType>& As = tSys->A();
    const meshField<Type,MeshType>& bs = tSys->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) = As(l,d,i,j,k)/s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) = bs(l,d,i,j,k)/s(l,d,i,j,k);

    if (shallow)
        collapse(s);

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const linearSystem<SType,Type,MeshType>& sys,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(sys/ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator/
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const tmp<meshField<scalar,MeshType>>& ts
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet(tSys/ts());

    if (ts.isTmp())
        ts.clear();

    return tRet;
}

// Summation

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator+
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        SysType::New
        (
            "(" + sys1.name() + "+" + sys2.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys1.x())
        );

    tRet.ref() = sys1;
    tRet.ref() += sys2;

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        reuseSysTmp<STypeR,SType1,Type,MeshType>::New
        (
            tSys1,
            "(" + tSys1->name() + "+" + sys2.name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != static_cast<const void*>(&tSys1()))
        tRet.ref() = tSys1;

    tRet.ref() += sys2;

    if (tSys1.isTmp())
        tSys1.clear();

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator+
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        reuseSysTmp<STypeR,SType1,Type,MeshType>::New
        (
            tSys2,
            "(" + sys1.name() + "+" + tSys2->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != static_cast<const void*>(&tSys2()))
        tRet.ref() = tSys2;
    tRet.ref() += sys1;

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        reuseSysTmpTmp<STypeR,SType1,SType1,SType2,Type,MeshType>::New
        (
            tSys1,
            tSys2,
            "(" + tSys1->name() + "+" + tSys2->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) == static_cast<const void*>(&tSys1()))
    {
        tRet.ref() += tSys2;
    }
    else if
    (
        static_cast<const void*>(&tRet())
     == static_cast<const void*>(&tSys2())
    )
    {
        tRet.ref() += tSys1;
    }
    else
    {
        if (std::is_same<STypeR,SType1>::value)
        {
            tRet.ref() = tSys1;
            tRet.ref() += tSys2;
        }
        else
        {
            tRet.ref() = tSys2;
            tRet.ref() += tSys1;
        }
    }

    if (tSys1.isTmp())
        tSys1.clear();

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const linearSystem<SType,Type,MeshType>& sys,
    const meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "+" + field.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() += field;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "+" + field.name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() += field;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const linearSystem<SType,Type,MeshType>& sys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "*" + tField->name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "+" + tField->name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() += tField();

    if (tSys.isTmp())
        tSys.clear();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const meshField<Type,MeshType>& field,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + field.name() + "+" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() += field;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const meshField<Type,MeshType>& field,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + field.name() + "+" + tSys->name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() += field;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<meshField<Type,MeshType>>& tField,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + tField->name() + "+" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<meshField<Type,MeshType>>& tField,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tField->name() + "+" + tSys->name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Subtraction

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator-
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        SysType::New
        (
            "(" + sys1.name() + "-" + sys2.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys1.x())
        );

    tRet.ref() = sys1;
    tRet.ref() -= sys2;

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet
    (
        reuseSysTmp<STypeR,SType1,Type,MeshType>::New
        (
            tSys1,
            "(" + tSys1->name() + "-" + sys2.name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != static_cast<const void*>(&tSys1()))
        tRet.ref() = tSys1;

    tRet.ref() -= sys2;

    if (tSys1.isTmp())
        tSys1.clear();

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator-
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet
    (
        reuseSysTmp<STypeR,SType1,Type,MeshType>::New
        (
            tSys2,
            "(" + sys1.name() + "-" + tSys2->name() + ")"
        )
    );

    if (static_cast<const void*>(&tRet()) != static_cast<const void*>(&tSys2()))
        tRet.ref() = tSys2;

    tRet.ref() *= -1;
    tRet.ref() += sys1;

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    typedef typename typeOfSum<SType1,SType2>::type STypeR;
    typedef linearSystem<STypeR,Type,MeshType> SysType;

    tmp<SysType> tRet =
        reuseSysTmpTmp<STypeR,SType1,SType1,SType2,Type,MeshType>::New
        (
            tSys1,
            tSys2,
            "(" + tSys1->name() + "-" + tSys2->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) == static_cast<const void*>(&tSys1()))
    {
        tRet.ref() -= tSys2;
    }
    else if
    (
        static_cast<const void*>(&tRet())
     == static_cast<const void*>(&tSys2())
    )
    {
        tRet.ref() *= -1;
        tRet.ref() += tSys1;
    }
    else
    {
        tRet.ref() = tSys1;
        tRet.ref() -= tSys2;
    }

    if (tSys1.isTmp())
        tSys1.clear();

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const linearSystem<SType,Type,MeshType>& sys,
    const meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "-" + field.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() -= field;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "-" + field.name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() -= field;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const linearSystem<SType,Type,MeshType>& sys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + sys.name() + "-" + tField->name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() -= tField();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tSys->name() + "-" + tField->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() -= tField();

    if (tSys.isTmp())
        tSys.clear();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const meshField<Type,MeshType>& field,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + field.name() + "-" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() *= -1;
    tRet.ref() += field;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const meshField<Type,MeshType>& field,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + field.name() + "-" + tSys->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() *= -1;
    tRet.ref() += field;

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<meshField<Type,MeshType>>& tField,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        linearSystem<SType,Type,MeshType>::New
        (
            "(" + tField->name() + "-" + sys.name() + ")",
            const_cast<meshField<Type,MeshType>&>(sys.x())
        );

    tRet.ref() = sys;
    tRet.ref() *= -1;
    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<meshField<Type,MeshType>>& tField,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet =
        reuseSysTmp<SType,SType,Type,MeshType>::New
        (
            tSys,
            "(" + tField->name() + "-" + tSys->name() + ")"
        );

    if (static_cast<const void*>(&tRet()) != &tSys())
        tRet.ref() = tSys;

    tRet.ref() *= -1;
    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Equation

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator==
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    return sys1 - sys2;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const linearSystem<SType2,Type,MeshType>& sys2
)
{
    return tSys1 - sys2;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator==
(
    const linearSystem<SType1,Type,MeshType>& sys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    return sys1 - tSys2;
}

template<class SType1, class SType2, class Type, class MeshType>
inline tmp<linearSystem<typename typeOfSum<SType1,SType2>::type,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys2
)
{
    return tSys1 - tSys2;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const linearSystem<SType,Type,MeshType>& sys,
    const meshField<Type,MeshType>& field
)
{
    return sys - field;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const meshField<Type,MeshType>& field
)
{
    return tSys - field;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const linearSystem<SType,Type,MeshType>& sys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    return sys - tField;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const tmp<meshField<Type,MeshType>>& tField
)
{
    return tSys - tField;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const meshField<Type,MeshType>& field,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    return field - sys;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const meshField<Type,MeshType>& field,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    return field - tSys;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<meshField<Type,MeshType>>& tField,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    return tField - sys;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<meshField<Type,MeshType>>& tField,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    return tField - tSys;
}

}

}

}
