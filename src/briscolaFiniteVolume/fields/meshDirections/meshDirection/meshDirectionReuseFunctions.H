#ifndef meshDirectionReuseFunctions_H
#define meshDirectionReuseFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template
<
    class TypeR,
    class Type1,
    class MeshType
>
class reuseDirectionTmp
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<Type1,MeshType>>& tf1
    )
    {
        return
            meshDirection<TypeR,MeshType>::New
            (
                tf1().fvMsh(),
                tf1().levelNum(),
                tf1().directionNum()
            );
    }
};

template
<
    class TypeR,
    class MeshType
>
class reuseDirectionTmp<TypeR, TypeR, MeshType>
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<TypeR,MeshType>>& tf1
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return
                meshDirection<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum(),
                    tf1->directionNum()
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
class reuseDirectionTmpTmp
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<Type1,MeshType>>& tf1,
        const tmp<meshDirection<Type2,MeshType>>& tf2
    )
    {
        return
            meshDirection<TypeR,MeshType>::New
            (
                tf1->fvMsh(),
                tf1->levelNum(),
                tf1->directionNum()
            );
    }
};

template
<
    class TypeR,
    class Type12,
    class Type2,
    class MeshType
>
class reuseDirectionTmpTmp<TypeR, TypeR, Type12, Type2, MeshType>
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<TypeR,MeshType>>& tf1,
        const tmp<meshDirection<Type2,MeshType>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return
                meshDirection<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum(),
                    tf1->directionNum()
                );
        }
    }
};

template
<
    class TypeR,
    class Type1,
    class Type12,
    class MeshType
>
class reuseDirectionTmpTmp<TypeR, Type1, Type12, TypeR, MeshType>
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<Type1,MeshType>>& tf1,
        const tmp<meshDirection<TypeR,MeshType>>& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return
                meshDirection<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum(),
                    tf1->directionNum()
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
class reuseDirectionTmpTmp<TypeR, TypeR, Type12, TypeR, MeshType>
{
public:

    static tmp<meshDirection<TypeR,MeshType>> New
    (
        const tmp<meshDirection<TypeR,MeshType>>& tf1,
        const tmp<meshDirection<TypeR,MeshType>>& tf2
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
            return
                meshDirection<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum(),
                    tf1->directionNum()
                );
        }
    }
};

}

}

}

#endif
