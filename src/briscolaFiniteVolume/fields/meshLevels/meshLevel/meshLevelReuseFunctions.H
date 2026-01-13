#ifndef meshLevelReuseFunctions_H
#define meshLevelReuseFunctions_H

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
class reuseLevelTmp
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<Type1,MeshType>>& tf1
    )
    {
        return
            meshLevel<TypeR,MeshType>::New
            (
                tf1->fvMsh(),
                tf1->levelNum()
            );
    }
};

template
<
    class TypeR,
    class MeshType
>
class reuseLevelTmp<TypeR, TypeR, MeshType>
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<TypeR,MeshType>>& tf1
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return
                meshLevel<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum()
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
class reuseLevelTmpTmp
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<Type1,MeshType>>& tf1,
        const tmp<meshLevel<Type2,MeshType>>& tf2
    )
    {
        return
            meshLevel<TypeR,MeshType>::New
            (
                tf1->fvMsh(),
                tf1->levelNum()
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
class reuseLevelTmpTmp<TypeR, TypeR, Type12, Type2, MeshType>
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<TypeR,MeshType>>& tf1,
        const tmp<meshLevel<Type2,MeshType>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return
                meshLevel<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum()
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
class reuseLevelTmpTmp<TypeR, Type1, Type12, TypeR, MeshType>
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<Type1,MeshType>>& tf1,
        const tmp<meshLevel<TypeR,MeshType>>& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return
                meshLevel<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum()
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
class reuseLevelTmpTmp<TypeR, TypeR, Type12, TypeR, MeshType>
{
public:

    static tmp<meshLevel<TypeR,MeshType>> New
    (
        const tmp<meshLevel<TypeR,MeshType>>& tf1,
        const tmp<meshLevel<TypeR,MeshType>>& tf2
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
                meshLevel<TypeR,MeshType>::New
                (
                    tf1->fvMsh(),
                    tf1->levelNum()
                );
        }
    }
};

}

}

}

#endif
