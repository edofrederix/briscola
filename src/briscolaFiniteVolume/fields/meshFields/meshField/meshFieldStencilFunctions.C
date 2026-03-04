#ifndef meshFieldStencilFunctions_C
#define meshFieldStencilFunctions_C

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void rowProduct
(
    meshField<typename stencilProduct<Type1,Type2>::type,MeshType>& res,
    const meshField<Type1,MeshType>& f1,
    const meshField<Type2,MeshType>& f2
)
{
    forAll(res, l)
        rowProduct(res[l], f1[l], f2[l]);
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshField<Type1,MeshType>& f1,
    const meshField<Type2,MeshType>& f2
)
{
    typedef typename stencilProduct<Type1,Type2>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + f1.name() + "," + f2.name() + ")",
            f1.fvMsh()
        );

    rowProduct(tRes.ref(), f1, f2);

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshField<Type1,MeshType>>& tf1,
    const meshField<Type2,MeshType>& f2
)
{
    if (tf1.isTmp())
        tf1->correctBoundaryConditions();

    typedef typename stencilProduct<Type1,Type2>::type productType;

    // Tmp cannot be reused otherwise aliasing may occur

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + tf1->name() + "," + f2.name() + ")",
            tf1->fvMsh()
        );

    rowProduct(tRes.ref(),tf1(),f2);

    if (tf1.isTmp())
        tf1.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshField<Type1,MeshType>& f1,
    const tmp<meshField<Type2,MeshType>>& tf2
)
{
    if (tf2.isTmp())
        tf2->correctBoundaryConditions();

    typedef typename stencilProduct<Type1,Type2>::type productType;

    // Tmp cannot be reused otherwise aliasing may occur

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + f1.name() + "," + tf2->name() + ")",
            tf2->fvMsh()
        );

    rowProduct(tRes.ref(),f1,tf2());

    if (tf2.isTmp())
        tf2.clear();

    return tRes;
}

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshField<Type1,MeshType>>& tf1,
    const tmp<meshField<Type2,MeshType>>& tf2
)
{
    if (tf1.isTmp())
        tf1->correctBoundaryConditions();

    if (tf2.isTmp())
        tf2->correctBoundaryConditions();

    typedef typename stencilProduct<Type1,Type2>::type productType;

    // Tmp cannot be reused otherwise aliasing may occur

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + tf1->name() + "," + tf2->name() + ")",
            tf1->fvMsh()
        );

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
    meshField<typename stencilProduct<Type,Form>::type,MeshType>& res,
    const meshField<Type,MeshType>& f1,
    const Form& s
)
{
    forAll(res, l)
        rowProduct(res[l], f1[l], s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const meshField<Type,MeshType>& f1, const Form& s)
{
    typedef typename stencilProduct<Type,Form>::type productType;

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + f1.name() + "," + Foam::name(s) + ")",
            f1.fvMsh()
        );

    rowProduct(tRes.ref(), f1, s);

    return tRes;
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const tmp<meshField<Type,MeshType>>& tf1, const Form& s)
{
    if (tf1.isTmp())
        tf1->correctBoundaryConditions();

    typedef typename stencilProduct<Type,Form>::type productType;

    // Tmp cannot be reused otherwise aliasing may occur

    tmp<meshField<productType,MeshType>> tRes =
        meshField<productType,MeshType>::New
        (
            "rowProduct(" + tf1->name() + "," + Foam::name(s) + ")",
            tf1->fvMsh()
        );

    rowProduct(tRes.ref(),tf1(),s);

    if (tf1.isTmp())
        tf1.clear();

    return tRes;
}

// Use rowProduct(s,A) = rowProduct(A,s)

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshField<typename stencilProduct<Form,Type>::type,MeshType>& res,
    const Form& s,
    const meshField<Type,MeshType>& f2
)
{
    rowProduct(res,f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form& s, const meshField<Type,MeshType>& f2)
{
    return rowProduct(f2,s);
}

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form& s, const tmp<meshField<Type,MeshType>>& tf2)
{
    return rowProduct(tf2,s);
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
    tmp<meshField<scalar,MeshType>> tRes =
        meshField<scalar,MeshType>::New
        (
            "rowSum(" + f.name() + ")",
            f.fvMsh()
        );

    rowSum(tRes.ref(), f);

    return tRes;
}

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
rowSum(const tmp<meshField<Type,MeshType>>& tf)
{
    if (tf.isTmp())
        tf->correctBoundaryConditions();

    tmp<meshField<scalar,MeshType>> tRes =
        meshField<scalar,MeshType>::New
        (
            "rowSum(" + tf->name() + ")",
            tf->fvMsh()
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
