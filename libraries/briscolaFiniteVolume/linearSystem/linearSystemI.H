namespace Foam
{

namespace briscola
{

namespace fv
{

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
    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

    tRet->A() *= s;
    tRet->b() *= s;

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
    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

    tRet->A() *= s;
    tRet->b() *= s;

    return tRet;
}

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
    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

    tRet->A() /= s;
    tRet->b() /= s;

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
operator*
(
    const meshField<scalar,MeshType>& s,
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    const_cast<meshField<scalar,MeshType>&>(s).restrict();

    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

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

    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

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

    tmp<linearSystem<SType,Type,MeshType>>& tRet =
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys);

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

}

}

}
