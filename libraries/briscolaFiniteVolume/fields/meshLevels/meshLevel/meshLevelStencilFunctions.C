#ifndef meshLevelStencilFunctions_C
#define meshLevelStencilFunctions_C

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>& res,
    const meshLevel<Type1,MeshType>& f1,
    const meshLevel<Type2,MeshType>& f2
)
{
    forAll(res, d)
        Amul(res[d], f1[d], f2[d]);
}

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshLevel<Type1,MeshType>& f1,
    const meshLevel<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes
    (
        new meshLevel<productType,MeshType>
        (
            f1.fvMsh(),
            f1.levelNum()
        )
    );

    Amul(tRes.ref(), f1, f2);

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshLevel<Type1,MeshType>>& tf1,
    const meshLevel<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes =
        reuseLevTmp<productType,Type1,MeshType>::New(tf1);

    Amul(tRes.ref(),tf1(),f2);

    tf1.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshLevel<Type1,MeshType>& f1,
    const tmp<meshLevel<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes =
        reuseLevTmp<productType,Type2,MeshType>::New(tf2);

    Amul(tRes.ref(),f1,tf2());

    tf2.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshLevel<Type1,MeshType>>& tf1,
    const tmp<meshLevel<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes =
        reuseLevTmpTmp<productType,Type1,Type1,Type2,MeshType>::New(tf1,tf2);

    Amul(tRes.ref(),tf1(),tf2());

    tf1.clear();
    tf2.clear();

    return tRes;
}

template<class Type, class Form, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Type,Form>::type,MeshType>& res,
    const meshLevel<Type,MeshType>& f1,
    const Form& s
)
{
    forAll(res, d)
        Amul(res[d], f1[d], s);
}

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshLevel<Type,MeshType>& f1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes
    (
        new meshLevel<productType,MeshType>
        (
            f1.fvMsh(),
            f1.levelNum()
        )
    );

    Amul(tRes.ref(), f1, s);

    return tRes;
}

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshLevel<Type,MeshType>>& tf1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshLevel<productType,MeshType>> tRes =
        reuseLevTmp<productType,Type,MeshType>::New(tf1);

    Amul(tRes.ref(),tf1(),s);

    tf1.clear();

    return tRes;
}

// Use Amul(s,A) = Amul(A,s)

template<class Type, class Form, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Form,Type>::type,MeshType>& res,
    const Form& s,
    const meshLevel<Type,MeshType>& f2
)
{
    Amul(res,f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const meshLevel<Type,MeshType>& f2)
{
    return Amul(f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const tmp<meshLevel<Type,MeshType>>& tf2)
{
    return Amul(tf2,s);
}

template<class Type, class MeshType>
void rowSum
(
    meshLevel<scalar,MeshType>& res,
    const meshLevel<Type,MeshType>& f
)
{
    forAll(res, d)
        rowSum(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
rowSum(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<scalar,MeshType>> tRes
    (
        new meshLevel<scalar,MeshType>
        (
            f.fvMsh(),
            f.levelNum()
        )
    );

    rowSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
rowSum(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<scalar,MeshType>> tRes
    (
        new meshLevel<scalar,MeshType>
        (
            tf->fvMsh(),
            tf->levelNum()
        )
    );

    rowSum(tRes.ref(), tf());

    tf.clear();

    return tRes;
}

template<class Type, class MeshType>
scalarList matrixSum(const meshLevel<Type,MeshType>& f)
{
    return gSum(rowSum(f));
}

template<class Type, class MeshType>
scalarList matrixSum(const tmp<meshLevel<Type,MeshType>>& tf)
{
    scalarList ret(matrixSum(tf()));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
void neighborSum
(
    meshLevel<scalar,MeshType>& res,
    const meshLevel<Type,MeshType>& f
)
{
    forAll(res, d)
        neighborSum(res[d], f[d]);
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
neighborSum(const meshLevel<Type,MeshType>& f)
{
    tmp<meshLevel<scalar,MeshType>> tRes
    (
        new meshLevel<scalar,MeshType>
        (
            f.fvMsh(),
            f.levelNum()
        )
    );

    neighborSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
neighborSum(const tmp<meshLevel<Type,MeshType>>& tf)
{
    tmp<meshLevel<scalar,MeshType>> tRes
    (
        new meshLevel<scalar,MeshType>
        (
            tf->fvMsh(),
            tf->levelNum()
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
