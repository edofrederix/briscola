namespace Foam
{

namespace briscola
{

namespace fv
{

// Multiplication

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator*
(
    const scalar s,
    const linearSystem<SType,Type,MeshType>& sys
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

    tRet->A() *= s;
    tRet->b() *= s;

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    tRet->A() *= s;
    tRet->b() *= s;

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

    tRet->A() *= s;
    tRet->b() *= s;

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    tRet->A() *= s;
    tRet->b() *= s;

    if (tSys.isTmp())
        tSys.clear();

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

    if (tSys.isTmp())
        tSys.clear();

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

    if (tSys.isTmp())
        tSys.clear();

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    tRet->A() /= s;
    tRet->b() /= s;

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) *= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) *= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) /= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) /= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

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
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    // The restricted field has undefined ghosts. Manually implement operation
    // to avoid segmentation errors.

    meshField<SType,MeshType>& A = tRet->A();
    meshField<Type,MeshType>& b = tRet->b();

    forAllCells(A, l, d, i, j, k)
        A(l,d,i,j,k) /= s(l,d,i,j,k);

    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) /= s(l,d,i,j,k);

    const_cast<meshField<scalar,MeshType>&>(s).makeShallow();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Summation

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys1)
    );

    tRet.ref() += sys2;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys1)
    );

    tRet.ref() += sys2;

    if (tSys1.isTmp())
        tSys1.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys2)
    );

    tRet.ref() += sys1;

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator+
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys1)
    );

    tRet.ref() += tSys2();

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Subtraction

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys1)
    );

    tRet.ref() -= sys2;

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys1)
    );

    tRet.ref() -= sys2;

    if (tSys1.isTmp())
        tSys1.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys2)
    );

    tRet.ref() *= -1;
    tRet.ref() += sys1;

    if (tSys2.isTmp())
        tSys2.clear();

    return tRet;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator-
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
)
{
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys1)
    );

    tRet.ref() -= tSys2();

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );


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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        new linearSystem<SType,Type,MeshType>(sys)
    );

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
    tmp<linearSystem<SType,Type,MeshType>> tRet
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys)
    );

    tRet.ref() *= -1;
    tRet.ref() += tField();

    if (tField.isTmp())
        tField.clear();

    if (tSys.isTmp())
        tSys.clear();

    return tRet;
}

// Equation

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    return sys1 - sys2;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const linearSystem<SType,Type,MeshType>& sys2
)
{
    return tSys1 - sys2;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const linearSystem<SType,Type,MeshType>& sys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
)
{
    return sys1 - tSys2;
}

template<class SType, class Type, class MeshType>
inline tmp<linearSystem<SType,Type,MeshType>>
operator==
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys1,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys2
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
