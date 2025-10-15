#ifndef meshFieldSetReuseFunctions_H
#define meshFieldSetReuseFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template
<
    class Type,
    class MeshType,
    direction N
>
class reuseSetTmp
{
public:

    static tmp<meshFieldSet<Type,MeshType,N>> New
    (
        const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
        const word name
    )
    {
        meshFieldSet<Type,MeshType,N>& set1 =
            const_cast<meshFieldSet<Type,MeshType,N>&>(tSet1());

        if (tSet1.isTmp())
        {
            set1.rename(name);
            return tSet1;
        }
        else
        {
            return
                meshFieldSet<Type,MeshType,N>::New
                (
                    name,
                    set1.fvMsh()
                );
        }
    }
};

template
<
    class Type,
    class MeshType,
    direction N
>
class reuseSetTmpTmp
{
public:

    static tmp<meshFieldSet<Type,MeshType,N>> New
    (
        const tmp<meshFieldSet<Type,MeshType,N>>& tSet1,
        const tmp<meshFieldSet<Type,MeshType,N>>& tSet2,
        const word name
    )
    {
        meshFieldSet<Type,MeshType,N>& set1 =
            const_cast<meshFieldSet<Type,MeshType,N>&>(tSet1());

        meshFieldSet<Type,MeshType,N>& set2 =
            const_cast<meshFieldSet<Type,MeshType,N>&>(tSet2());

        if (tSet1.isTmp())
        {
            set1.rename(name);
            return tSet1;
        }
        else if (tSet2.isTmp())
        {
            set2.rename(name);
            return tSet2;
        }
        else
        {
            return
                meshFieldSet<Type,MeshType,N>::New
                (
                    name,
                    set1.fvMsh()
                );
        }
    }
};

}

}

}

#endif
