#ifndef meshDirectionStencilFunctions_C
#define meshDirectionStencilFunctions_C

#include "meshDirectionReuseFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>& res,
    const meshDirection<Type1,MeshType>& f1,
    const meshDirection<Type2,MeshType>& f2
)
{
    res = Zero;

    // Make sure we iterate over the field and not the stencil, because the
    // field may have inactive cells due to boundary conditions.

    if (isStencil<Type1>())
    {
        forAllCells(f2, i, j, k)
            res(i,j,k) = Amul(f1,f2,i,j,k);
    }
    else
    {
        forAllCells(f1, i, j, k)
            res(i,j,k) = Amul(f1,f2,i,j,k);
    }
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshDirection<Type1,MeshType>& f1,
    const meshDirection<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes
    (
        new meshDirection<productType,MeshType>
        (
            f1.fvMsh(),
            f1.levelNum(),
            f1.directionNum()
        )
    );

    Amul(tRes.ref(), f1, f2);

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshDirection<Type1,MeshType>>& tf1,
    const meshDirection<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirTmp<productType,Type1,MeshType>::New(tf1);

    Amul(tRes.ref(),tf1(),f2);

    tf1.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshDirection<Type1,MeshType>& f1,
    const tmp<meshDirection<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirTmp<productType,Type2,MeshType>::New(tf2);

    Amul(tRes.ref(),f1,tf2());

    tf2.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshDirection<Type1,MeshType>>& tf1,
    const tmp<meshDirection<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirTmpTmp<productType,Type1,Type1,Type2,MeshType>::New(tf1,tf2);

    Amul(tRes.ref(),tf1(),tf2());

    tf1.clear();
    tf2.clear();

    return tRes;
}

template<class Type, class Form, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Type,Form>::type,MeshType>& res,
    const meshDirection<Type,MeshType>& f1,
    const Form& s
)
{
    res = Zero;

    forAllCells(f1, i, j, k)
    {
        res(i,j,k) = Amul(f1,s,i,j,k);
    }
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshDirection<Type,MeshType>& f1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes
    (
        new meshDirection<productType,MeshType>
        (
            f1.fvMsh(),
            f1.levelNum(),
            f1.directionNum()
        )
    );

    Amul(tRes.ref(), f1, s);

    return tRes;
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshDirection<Type,MeshType>>& tf1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirTmp<productType,Type,MeshType>::New(tf1);

    Amul(tRes.ref(),tf1(),s);

    tf1.clear();

    return tRes;
}

// Use Amul(s,A) = Amul(A,s)

template<class Type, class Form, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Form,Type>::type,MeshType>& res,
    const Form& s,
    const meshDirection<Type,MeshType>& f2
)
{
    Amul(res, f2, s);
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const meshDirection<Type,MeshType>& f2)
{
    return Amul(f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const tmp<meshDirection<Type,MeshType>>& tf2)
{
    return Amul(tf2,s);
}

template<class Type, class MeshType>
void rowSum
(
    meshDirection<scalar,MeshType>& res,
    const meshDirection<Type,MeshType>& f
)
{
    res = Zero;

    forAllCells(f, i, j, k)
        res(i,j,k) = stencilSum(f(i,j,k));
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
rowSum(const meshDirection<Type,MeshType>& f)
{
    tmp<meshDirection<scalar,MeshType>> tRes
    (
        new meshDirection<scalar,MeshType>
        (
            f.fvMsh(),
            f.levelNum(),
            f.directionNum()
        )
    );

    rowSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
rowSum(const tmp<meshDirection<Type,MeshType>>& tf)
{
    tmp<meshDirection<scalar,MeshType>> tRes
    (
        new meshDirection<scalar,MeshType>
        (
            tf->fvMsh(),
            tf->levelNum(),
            tf->directionNum()
        )
    );

    rowSum(tRes.ref(), tf());

    tf.clear();

    return tRes;
}

template<class Type, class MeshType>
scalar matrixSum(const meshDirection<Type,MeshType>& f)
{
    return gSum(rowSum(f));
}

template<class Type, class MeshType>
scalar matrixSum(const tmp<meshDirection<Type,MeshType>>& tf)
{
    scalar ret(matrixSum(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
void neighborSum
(
    meshDirection<scalar,MeshType>& res,
    const meshDirection<Type,MeshType>& f
)
{
    res = Zero;

    forAllCells(f, i, j, k)
        res(i,j,k) = neighborSum(f(i,j,k));
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
neighborSum(const meshDirection<Type,MeshType>& f)
{
    tmp<meshDirection<scalar,MeshType>> tRes
    (
        new meshDirection<scalar,MeshType>
        (
            f.fvMsh(),
            f.levelNum(),
            f.directionNum()
        )
    );

    neighborSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
neighborSum(const tmp<meshDirection<Type,MeshType>>& tf)
{
    tmp<meshDirection<scalar,MeshType>> tRes
    (
        new meshDirection<scalar,MeshType>
        (
            tf->fvMsh(),
            tf->levelNum(),
            tf->directionNum()
        )
    );

    neighborSum(tRes.ref(), tf());

    tf.clear();

    return tRes;
}

}

}

}

#endif
