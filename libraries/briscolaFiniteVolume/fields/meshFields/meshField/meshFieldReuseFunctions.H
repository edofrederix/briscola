#ifndef meshFieldReuseFunctions_H
#define meshFieldReuseFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template
<
    class TypeR,
    class MeshType
>
tmp<meshField<TypeR,MeshType>> New
(
    const tmp<meshField<TypeR,MeshType>>& tf1,
    const word name,
    const bool initRet = false
)
{
    meshField<TypeR,MeshType>& f1 =
        const_cast<meshField<TypeR,MeshType>&>(tf1());

    if (tf1.isTmp())
    {
        f1.rename(name);
        return tf1;
    }
    else
    {
        tmp<meshField<TypeR,MeshType>> rtf
        (
            new meshField<TypeR,MeshType>(name, f1.fvMsh())
        );

        if (initRet)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
}

template
<
    class TypeR,
    class Type1,
    class MeshType
>
class reuseFieTmp
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<Type1,MeshType>>& tf1,
        const word name
    )
    {
        return tmp<meshField<TypeR,MeshType>>
        (
            new meshField<TypeR,MeshType>(name, tf1->fvMsh())
        );
    }
};

template
<
    class TypeR,
    class MeshType
>
class reuseFieTmp<TypeR, TypeR, MeshType>
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<TypeR,MeshType>>& tf1,
        const word name
    )
    {
        meshField<TypeR,MeshType>& f1 =
            const_cast<meshField<TypeR,MeshType>&>(tf1());

        if (tf1.isTmp())
        {
            f1.rename(name);
            return tf1;
        }
        else
        {
            return tmp<meshField<TypeR,MeshType>>
            (
                new meshField<TypeR,MeshType>(name, f1.fvMsh())
            );
        }
    }
};

template
<
    class TypeR,
    class Type1,
    class Type12,
    class Type2,
    class MeshType
>
class reuseFieTmpTmp
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<Type1,MeshType>>& tf1,
        const tmp<meshField<Type2,MeshType>>& tf2,
        const word name
    )
    {
        return tmp<meshField<TypeR,MeshType>>
        (
            new meshField<TypeR,MeshType>(name, tf1->fvMsh())
        );
    }
};

template
<
    class TypeR,
    class Type1,
    class Type12,
    class MeshType
>
class reuseFieTmpTmp<TypeR, Type1, Type12, TypeR, MeshType>
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<Type1,MeshType>>& tf1,
        const tmp<meshField<TypeR,MeshType>>& tf2,
        const word name
    )
    {
        meshField<TypeR,MeshType>& f2 =
            const_cast<meshField<TypeR,MeshType>&>(tf2());

        if (tf2.isTmp())
        {
            f2.rename(name);
            return tf2;
        }
        else
        {
            return tmp<meshField<TypeR,MeshType>>
            (
                new meshField<TypeR,MeshType>(name, f2.fvMsh())
            );
        }
    }
};

template
<
    class TypeR,
    class Type12,
    class Type2,
    class MeshType
>
class reuseFieTmpTmp<TypeR, TypeR, Type12, Type2, MeshType>
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<TypeR,MeshType>>& tf1,
        const tmp<meshField<Type2,MeshType>>& tf2,
        const word name
    )
    {
        meshField<TypeR,MeshType>& f1 =
            const_cast<meshField<TypeR,MeshType>&>(tf1());

        if (tf1.isTmp())
        {
            f1.rename(name);
            return tf1;
        }
        else
        {
            return tmp<meshField<TypeR,MeshType>>
            (
                new meshField<TypeR,MeshType>(name, f1.fvMsh())
            );
        }
    }
};

template
<
    class TypeR,
    class Type12,
    class MeshType
>
class reuseFieTmpTmp<TypeR, TypeR, Type12, TypeR, MeshType>
{
public:

    static tmp<meshField<TypeR,MeshType>> New
    (
        const tmp<meshField<TypeR,MeshType>>& tf1,
        const tmp<meshField<TypeR,MeshType>>& tf2,
        const word name
    )
    {
        meshField<TypeR,MeshType>& f1 =
            const_cast<meshField<TypeR,MeshType>&>(tf1());

        meshField<TypeR,MeshType>& f2 =
            const_cast<meshField<TypeR,MeshType>&>(tf2());

        if (tf1.isTmp())
        {
            f1.rename(name);
            return tf1;
        }
        else if (tf2.isTmp())
        {
            f2.rename(name);
            return tf2;
        }
        else
        {
            return tmp<meshField<TypeR,MeshType>>
            (
                new meshField<TypeR,MeshType>(name, f1.fvMsh())
            );
        }
    }
};

}

}

}

#endif
