#ifndef meshFieldStencilFunctions_C
#define meshFieldStencilFunctions_C

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void Amul
(
    meshField<typename stencilProduct<Type1,Type2>::type,MeshType>& res,
    const meshField<Type1,MeshType>& f1,
    const meshField<Type2,MeshType>& f2
)
{
    forAll(res, l)
        Amul(res[l], f1[l], f2[l]);
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshField<Type1,MeshType>& f1,
    const meshField<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshField<productType,MeshType>> tRes
    (
        new meshField<productType,MeshType>
        (
            "Amul(" + f1.name() + "," + f2.name() + ")",
            f1.fvMsh()
        )
    );

    Amul(tRes.ref(), f1, f2);

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshField<Type1,MeshType>>& tf1,
    const meshField<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        reuseFieTmp<productType,Type1,MeshType>::New
        (
            tf1,
            "Amul(" + tf1->name() + "," + f2.name() + ")"
        );

    Amul(tRes.ref(),tf1(),f2);

    tf1.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshField<Type1,MeshType>& f1,
    const tmp<meshField<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        reuseFieTmp<productType,Type2,MeshType>::New
        (
            tf2,
            "Amul(" + f1.name() + "," + tf2->name() + ")"
        );

    Amul(tRes.ref(),f1,tf2());

    tf2.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshField<Type1,MeshType>>& tf1,
    const tmp<meshField<Type2,MeshType>>& tf2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        reuseFieTmpTmp<productType,Type1,Type1,Type2,MeshType>::New
        (
            tf1,
            tf2,
            "Amul(" + tf1->name() + "," + tf2->name() + ")"
        );

    Amul(tRes.ref(),tf1(),tf2());

    tf1.clear();
    tf2.clear();

    return tRes;
}

template<class Type, class Form, class MeshType>
void Amul
(
    meshField<typename stencilProduct<Type,Form>::type,MeshType>& res,
    const meshField<Type,MeshType>& f1,
    const Form& s
)
{
    forAll(res, l)
        Amul(res[l], f1[l], s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshField<Type,MeshType>& f1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshField<productType,MeshType>> tRes
    (
        new meshField<productType,MeshType>
        (
            "Amul(" + f1.name() + "," + Foam::name(s) + ")",
            f1.fvMsh()
        )
    );

    Amul(tRes.ref(), f1, s);

    return tRes;
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshField<Type,MeshType>>& tf1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        reuseFieTmp<productType,Type,MeshType>::New
        (
            tf1,
            "Amul(" + tf1->name() + "," + Foam::name(s) + ")"
        );

    Amul(tRes.ref(),tf1(),s);

    tf1.clear();

    return tRes;
}

// Use Amul(s,A) = Amul(A,s)

template<class Type, class Form, class MeshType>
void Amul
(
    meshField<typename stencilProduct<Form,Type>::type,MeshType>& res,
    const Form& s,
    const meshField<Type,MeshType>& f2
)
{
    Amul(res,f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const meshField<Type,MeshType>& f2)
{
    return Amul(f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form& s, const tmp<meshField<Type,MeshType>>& tf2)
{
    return Amul(tf2,s);
}

template<class Type, class MeshType>
void rowSum
(
    meshField<scalar,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, l)
        rowSum(res[l], f[l]);
}

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
rowSum(const meshField<Type,MeshType>& f)
{
    tmp<meshField<scalar,MeshType>> tRes
    (
        new meshField<scalar,MeshType>
        (
            "rowSum(" + f.name() + ")",
            f.fvMsh()
        )
    );

    rowSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
rowSum(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<scalar,MeshType>> tRes
    (
        new meshField<scalar,MeshType>
        (
            "rowSum(" + tf->name() + ")",
            tf->fvMsh()
        )
    );

    rowSum(tRes.ref(), tf());

    tf.clear();

    return tRes;
}

template<class Type, class MeshType>
scalarList matrixSum(const meshField<Type,MeshType>& f)
{
    return gSum(rowSum(f[0]));
}

template<class Type, class MeshType>
scalarList matrixSum(const tmp<meshField<Type,MeshType>>& tf)
{
    scalarList ret(matrixSum(tf()[0]));
    tf.clear();
    return ret;
}

template<class Type, class MeshType>
void neighborSum
(
    meshField<scalar,MeshType>& res,
    const meshField<Type,MeshType>& f
)
{
    forAll(res, l)
        neighborSum(res[l], f[l]);
}

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
neighborSum(const meshField<Type,MeshType>& f)
{
    tmp<meshField<scalar,MeshType>> tRes
    (
        new meshField<scalar,MeshType>
        (
            "neighborSum(" + f.name() + ")",
            f.fvMsh()
        )
    );

    neighborSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
neighborSum(const tmp<meshField<Type,MeshType>>& tf)
{
    tmp<meshField<scalar,MeshType>> tRes
    (
        new meshField<scalar,MeshType>
        (
            "neighborSum(" + tf->name() + ")",
            tf->fvMsh()
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
