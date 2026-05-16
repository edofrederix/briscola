#include "PstreamReduceOps.H"
#include "meshDirectionReuseFunctions.H"

#define TEMPLATE template<class Type, class MeshType>
#include "meshDirectionFunctionsM.C"

// Scalar return type must be deduced because of cell space
#define SCALARPRODTYPE typename scalarProduct<Type,Type>::type

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void mag
(
    meshDirection<SCALARPRODTYPE,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    mag(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<SCALARPRODTYPE,MeshType>>
mag(const meshDirection<Type,MeshType>& D)
{
    tmp<meshDirection<SCALARPRODTYPE,MeshType>> tRes =
        meshDirection<SCALARPRODTYPE,MeshType>::New
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        );

    mag(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<SCALARPRODTYPE,MeshType>>
mag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    tmp<meshDirection<SCALARPRODTYPE,MeshType>> tRes =
        reuseDirectionTmp<SCALARPRODTYPE,Type,MeshType>::New(tD);

    mag(tRes.ref(), tD());
    if (tD.isTmp())
        tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptMag
(
    meshDirection<Type,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptMag(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const meshDirection<Type,MeshType>& D)
{
    tmp<meshDirection<Type,MeshType>> tRes =
        meshDirection<Type,MeshType>::New
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        );

    cmptMag(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptMag(const tmp<meshDirection<Type,MeshType>>& tD)
{
    tmp<meshDirection<Type,MeshType>> tRes =
        reuseDirectionTmp<Type,Type,MeshType>::New(tD);

    cmptMag(tRes.ref(), tD());
    if (tD.isTmp())
        tD.clear();
    return tRes;
}

template<class Type, class MeshType>
void cmptSqr
(
    meshDirection<Type,MeshType>& res,
    const meshDirection<Type,MeshType>& D
)
{
    cmptSqr(res.B(), D.B());
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptSqr(const meshDirection<Type,MeshType>& D)
{
    tmp<meshDirection<Type,MeshType>> tRes =
        meshDirection<Type,MeshType>::New
        (
            D.fvMsh(),
            D.levelNum(),
            D.directionNum()
        );

    cmptSqr(tRes.ref(), D);
    return tRes;
}

template<class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
cmptSqr(const tmp<meshDirection<Type,MeshType>>& tD)
{
    tmp<meshDirection<Type,MeshType>> tRes =
        reuseDirectionTmp<Type,Type,MeshType>::New(tD);

    cmptSqr(tRes.ref(), tD());
    if (tD.isTmp())
        tD.clear();
    return tRes;
}

template<class Type, class MeshType>
Type max(const meshDirection<Type,MeshType>& D)
{
    Type Max(pTraits<Type>::min);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (D(i,j,k) > Max)
        {
            Max = D(i,j,k);
        }
    }

    return Max;
}

template<class Type, class MeshType>
Type max(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(max(tD()));
    if (tD.isTmp())
        tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type min(const meshDirection<Type,MeshType>& D)
{
    Type Min(pTraits<Type>::max);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        if (D(i,j,k) < Min)
        {
            Min = D(i,j,k);
        }
    }

    return Min;
}

template<class Type, class MeshType>
Type min(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(min(tD()));
    if (tD.isTmp())
        tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type sum(const meshDirection<Type,MeshType>& D)
{
    Type Sum(Zero);

    // Only for interior cells

    forAllCells(D, i, j, k)
    {
        Sum += D(i,j,k);
    }

    return Sum;
}

template<class Type, class MeshType>
Type sum(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(sum(tD()));
    if (tD.isTmp())
        tD.clear();
    return ret;
}

template<class Type, class MeshType>
Type average(const meshDirection<Type,MeshType>& D)
{
    return sum(D)/D.size();
}

template<class Type, class MeshType>
Type average(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(average(tD()));
    if (tD.isTmp())
        tD.clear();
    return ret;
}

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type, class MeshType>                                           \
ReturnType gFunc(const meshDirection<Type,MeshType>& D)                        \
{                                                                              \
    ReturnType res(Func(D));                                                   \
    const label comm = D.lvl().comms();                                        \
    reduce(res, rFunc##Op<ReturnType>(), Pstream::msgType(), comm);            \
    return res;                                                                \
}                                                                              \
                                                                               \
template<class Type, class MeshType>                                           \
ReturnType gFunc(const tmp<meshDirection<Type,MeshType>>& tD)                  \
{                                                                              \
    ReturnType ret(gFunc(tD()));                                               \
    if (tD.isTmp()) tD.clear();                                                \
    return ret;                                                                \
}

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)

#undef G_UNARY_FUNCTION

template<class Type, class MeshType>
Type gAverage(const meshDirection<Type,MeshType>& D)
{
    label n = D.size();
    Type s = sum(D);
    sumReduce(s, n, Pstream::msgType(), D.lvl().comms());

    if (n > 0)
    {
        return s/n;
    }
    else
    {
        return Zero;
    }
}

template<class Type, class MeshType>
Type gAverage(const tmp<meshDirection<Type,MeshType>>& tD)
{
    Type ret(gAverage(tD()));
    if (tD.isTmp())
        tD.clear();
    return ret;
}

BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<typename product<Type1, Type2>::type,MeshType>& res,         \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D1.B(), D2.B());                                           \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>              \
operator Op                                                                    \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        meshDirection<productType,MeshType>::New                               \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), D1, D2);                                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>              \
operator Op                                                                    \
(                                                                              \
    const meshDirection<Type1,MeshType>& D1,                                   \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type2,MeshType>::New(tD2);               \
    OpFunc(tRes.ref(), D1, tD2());                                             \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>              \
operator Op                                                                    \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const meshDirection<Type2,MeshType>& D2                                    \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type1,MeshType>::New(tD1);               \
    OpFunc(tRes.ref(), tD1(), D2);                                             \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2, class MeshType>                             \
tmp<meshDirection<typename product<Type1, Type2>::type,MeshType>>              \
operator Op                                                                    \
(                                                                              \
    const tmp<meshDirection<Type1,MeshType>>& tD1,                             \
    const tmp<meshDirection<Type2,MeshType>>& tD2                              \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmpTmp<productType,Type1,Type1,Type2,MeshType>::         \
        New(tD1, tD2);                                                         \
    OpFunc(tRes.ref(), tD1(), tD2());                                          \
    if (tD1.isTmp()) tD1.clear();                                              \
    if (tD2.isTmp()) tD2.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
/* VectorSpace */                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<typename product<Type, Form>::type,MeshType>& res,           \
    const meshDirection<Type,MeshType>& D1,                                    \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D1.B(), vs);                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const meshDirection<Type,MeshType>& D1,                                    \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        meshDirection<productType,MeshType>::New                               \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), D1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const tmp<meshDirection<Type,MeshType>>& tD1,                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), tD1(), static_cast<const Form&>(vs));                   \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<typename product<Form, Type>::type,MeshType>& res,           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const meshDirection<Type,MeshType>& D1                                     \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), vs, D1.B());                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const meshDirection<Type,MeshType>& D1                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        meshDirection<productType,MeshType>::New                               \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), D1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<meshDirection<Type,MeshType>>& tD1                               \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tD1());                   \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
/* CellSpace */                                                                \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<typename product<Type, Form>::type,MeshType>& res,           \
    const meshDirection<Type,MeshType>& D1,                                    \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), D1.B(), vs);                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const meshDirection<Type,MeshType>& D1,                                    \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        meshDirection<productType,MeshType>::New                               \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), D1, static_cast<const Form&>(vs));                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Type, Form>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const tmp<meshDirection<Type,MeshType>>& tD1,                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), tD1(), static_cast<const Form&>(vs));                   \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
void OpFunc                                                                    \
(                                                                              \
    meshDirection<typename product<Form, Type>::type,MeshType>& res,           \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const meshDirection<Type,MeshType>& D1                                     \
)                                                                              \
{                                                                              \
    OpFunc(res.B(), vs, D1.B());                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const meshDirection<Type,MeshType>& D1                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        meshDirection<productType,MeshType>::New                               \
        (                                                                      \
            D1.fvMsh(),                                                        \
            D1.levelNum(),                                                     \
            D1.directionNum()                                                  \
        );                                                                     \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), D1);                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Type,                                                                \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class MeshType                                                             \
>                                                                              \
tmp<meshDirection<typename product<Form, Type>::type,MeshType>>                \
operator Op                                                                    \
(                                                                              \
    const CellSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const tmp<meshDirection<Type,MeshType>>& tD1                               \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    tmp<meshDirection<productType,MeshType>> tRes =                            \
        reuseDirectionTmp<productType,Type,MeshType>::New(tD1);                \
    OpFunc(tRes.ref(), static_cast<const Form&>(vs), tD1());                   \
    if (tD1.isTmp()) tD1.clear();                                              \
    return tRes;                                                               \
}

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR

}

}

}

#undef SCALARPRODTYPE

#include "undefBlockFunctionsM.H"
