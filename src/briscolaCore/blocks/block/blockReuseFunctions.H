#ifndef blockReuseFunctions_H
#define blockReuseFunctions_H

namespace Foam
{

namespace briscola
{

template<class TypeR, int P>
tmp<block<TypeR,P>> New
(
    const tmp<block<TypeR,P>>& tf1,
    const bool initRet = false
)
{
    if (tf1.isTmp())
    {
        return tf1;
    }
    else
    {
        tmp<block<TypeR,P>> rtf(new block<TypeR,P>(tf1().shape()));

        if (initRet)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
}

template<class TypeR, class Type1, int P>
class reuseTmp
{
public:

    static tmp<block<TypeR,P>> New(const tmp<block<Type1,P>>& tf1)
    {
        return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
    }
};

template<class TypeR, int P>
class reuseTmp<TypeR, TypeR, P>
{
public:

    static tmp<block<TypeR,P>> New(const tmp<block<TypeR,P>>& tf1)
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
        }
    }
};

template<class TypeR, class Type1, class Type12, class Type2, int P>
class reuseTmpTmp
{
public:

    static tmp<block<TypeR,P>> New
    (
        const tmp<block<Type1,P>>& tf1,
        const tmp<block<Type2,P>>& tf2
    )
    {
        return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
    }
};

template<class TypeR, class Type1, class Type12, int P>
class reuseTmpTmp<TypeR, Type1, Type12, TypeR, P>
{
public:

    static tmp<block<TypeR,P>> New
    (
        const tmp<block<Type1,P>>& tf1,
        const tmp<block<TypeR,P>>& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
        }
    }
};

template<class TypeR, class Type2, int P>
class reuseTmpTmp<TypeR, TypeR, TypeR, Type2, P>
{
public:

    static tmp<block<TypeR,P>> New
    (
        const tmp<block<TypeR,P>>& tf1,
        const tmp<block<Type2,P>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
        }
    }
};

template<class TypeR, int P>
class reuseTmpTmp<TypeR, TypeR, TypeR, TypeR, P>
{
public:

    static tmp<block<TypeR,P>> New
    (
        const tmp<block<TypeR,P>>& tf1,
        const tmp<block<TypeR,P>>& tf2
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
            return tmp<block<TypeR,P>>(new block<TypeR,P>(tf1().shape()));
        }
    }
};

}

}

#endif
