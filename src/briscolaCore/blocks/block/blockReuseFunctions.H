#ifndef blockReuseFunctions_H
#define blockReuseFunctions_H

namespace Foam
{

namespace briscola
{

template<class TypeR>
tmp<block<TypeR>> New
(
    const tmp<block<TypeR>>& tf1,
    const bool initRet = false
)
{
    if (tf1.isTmp())
    {
        return tf1;
    }
    else
    {
        tmp<block<TypeR>> rtf(new block<TypeR>(tf1().shape()));

        if (initRet)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
}

template<class TypeR, class Type1>
class reuseTmp
{
public:

    static tmp<block<TypeR>> New(const tmp<block<Type1>>& tf1)
    {
        return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
    }
};

template<class TypeR>
class reuseTmp<TypeR, TypeR>
{
public:

    static tmp<block<TypeR>> New(const tmp<block<TypeR>>& tf1)
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
        }
    }
};

template<class TypeR, class Type1, class Type12, class Type2>
class reuseTmpTmp
{
public:

    static tmp<block<TypeR>> New
    (
        const tmp<block<Type1>>& tf1,
        const tmp<block<Type2>>& tf2
    )
    {
        return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
    }
};

template<class TypeR, class Type1, class Type12>
class reuseTmpTmp<TypeR, Type1, Type12, TypeR>
{
public:

    static tmp<block<TypeR>> New
    (
        const tmp<block<Type1>>& tf1,
        const tmp<block<TypeR>>& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
        }
    }
};

template<class TypeR, class Type2>
class reuseTmpTmp<TypeR, TypeR, TypeR, Type2>
{
public:

    static tmp<block<TypeR>> New
    (
        const tmp<block<TypeR>>& tf1,
        const tmp<block<Type2>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
        }
    }
};

template<class TypeR>
class reuseTmpTmp<TypeR, TypeR, TypeR, TypeR>
{
public:

    static tmp<block<TypeR>> New
    (
        const tmp<block<TypeR>>& tf1,
        const tmp<block<TypeR>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return tmp<block<TypeR>>(new block<TypeR>(tf1().shape()));
        }
    }
};

}

}

#endif
