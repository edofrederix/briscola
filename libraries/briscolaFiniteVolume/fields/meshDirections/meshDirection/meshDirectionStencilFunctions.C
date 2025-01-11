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
void rowProduct
(
    meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>& res,
    const meshDirection<Type1,MeshType>& f1,
    const meshDirection<Type2,MeshType>& f2
)
{
    res = Zero;

    if (isStencil<Type1>())
    {
        forAllCells(f2, i, j, k)
            res(i,j,k) = rowProduct(f1,f2,i,j,k);
    }
    else
    {
        forAllCells(f1, i, j, k)
            res(i,j,k) = rowProduct(f1,f2,i,j,k);
    }
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
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

    rowProduct(tRes.ref(), f1, f2);

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshDirection<Type1,MeshType>>& tf1,
    const meshDirection<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirectionTmp<productType,Type1,MeshType>::New(tf1);

    rowProduct(tRes.ref(),tf1(),f2);

    if (tf1.isTmp())
        tf1.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshDirection<Type1,MeshType>& f1,
    const tmp<meshDirection<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirectionTmp<productType,Type2,MeshType>::New(tf2);

    rowProduct(tRes.ref(),f1,tf2());

    if (tf2.isTmp())
        tf2.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshDirection<Type1,MeshType>>& tf1,
    const tmp<meshDirection<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirectionTmpTmp<productType,Type1,Type1,Type2,MeshType>::New(tf1,tf2);

    rowProduct(tRes.ref(),tf1(),tf2());

    if (tf1.isTmp())
        tf1.clear();
    if (tf2.isTmp())
        tf2.clear();

    return tRes;
}

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshDirection<typename stencilProduct<Type,Form>::type,MeshType>& res,
    const meshDirection<Type,MeshType>& f1,
    const Form& s
)
{
    res = Zero;

    forAllCells(f1, i, j, k)
        res(i,j,k) = rowProduct(f1,s,i,j,k);
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const meshDirection<Type,MeshType>& f1, const Form& s)
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

    rowProduct(tRes.ref(), f1, s);

    return tRes;
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const tmp<meshDirection<Type,MeshType>>& tf1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshDirection<productType,MeshType>> tRes =
        reuseDirectionTmp<productType,Type,MeshType>::New(tf1);

    rowProduct(tRes.ref(),tf1(),s);

    if (tf1.isTmp())
        tf1.clear();

    return tRes;
}

// Use rowProduct(s,A) = rowProduct(A,s)

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshDirection<typename stencilProduct<Form,Type>::type,MeshType>& res,
    const Form& s,
    const meshDirection<Type,MeshType>& f2
)
{
    rowProduct(res, f2, s);
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form& s, const meshDirection<Type,MeshType>& f2)
{
    return rowProduct(f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form& s, const tmp<meshDirection<Type,MeshType>>& tf2)
{
    return rowProduct(tf2,s);
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
        res(i,j,k) = rowSum(f,i,j,k);
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

    if (tf.isTmp())
        tf.clear();

    return tRes;
}

}

}

}

#endif
