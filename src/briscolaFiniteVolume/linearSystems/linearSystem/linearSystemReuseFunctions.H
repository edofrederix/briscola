#ifndef linearSystemReuseFunctions_H
#define linearSystemReuseFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template
<
    class STypeR,
    class SType1,
    class Type,
    class MeshType
>
class reuseSysTmp
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
        const word name
    )
    {
        return tmp<linearSystem<STypeR,Type,MeshType>>
        (
            linearSystem<STypeR,Type,MeshType>::New
            (
                name,
                tSys1->x()
            )
        );
    }
};

template
<
    class STypeR,
    class Type,
    class MeshType
>
class reuseSysTmp<STypeR, STypeR, Type, MeshType>
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<STypeR,Type,MeshType>>& tSys1,
        const word name
    )
    {
        linearSystem<STypeR,Type,MeshType>& sys1 =
            const_cast<linearSystem<STypeR,Type,MeshType>&>(tSys1());

        if (tSys1.isTmp())
        {
            sys1.rename(name);
            return tSys1;
        }
        else
        {
            return tmp<linearSystem<STypeR,Type,MeshType>>
            (
                linearSystem<STypeR,Type,MeshType>::New
                (
                    name,
                    sys1.x()
                )
            );
        }
    }
};

template
<
    class STypeR,
    class SType1,
    class SType12,
    class SType2,
    class Type,
    class MeshType
>
class reuseSysTmpTmp
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
        const tmp<linearSystem<SType2,Type,MeshType>>& tSys2,
        const word name
    )
    {
        return tmp<linearSystem<STypeR,Type,MeshType>>
        (
            linearSystem<STypeR,Type,MeshType>::New
            (
                name,
                tSys1->x()
            )
        );
    }
};

template
<
    class STypeR,
    class SType1,
    class SType12,
    class Type,
    class MeshType
>
class reuseSysTmpTmp<STypeR, SType1, SType12, STypeR, Type, MeshType>
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<SType1,Type,MeshType>>& tSys1,
        const tmp<linearSystem<STypeR,Type,MeshType>>& tSys2,
        const word name
    )
    {
        linearSystem<STypeR,Type,MeshType>& sys2 =
            const_cast<linearSystem<STypeR,Type,MeshType>&>(tSys2());

        if (tSys2.isTmp())
        {
            sys2.rename(name);
            return tSys2;
        }
        else
        {
            return tmp<linearSystem<STypeR,Type,MeshType>>
            (
                linearSystem<STypeR,Type,MeshType>::New
                (
                    name,
                    sys2.x()
                )
            );
        }
    }
};

template
<
    class STypeR,
    class SType12,
    class SType2,
    class Type,
    class MeshType
>
class reuseSysTmpTmp<STypeR, STypeR, SType12, SType2, Type, MeshType>
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<STypeR,Type,MeshType>>& tSys1,
        const tmp<linearSystem<SType2,Type,MeshType>>& tSys2,
        const word name
    )
    {
        linearSystem<STypeR,Type,MeshType>& sys1 =
            const_cast<linearSystem<STypeR,Type,MeshType>&>(tSys1());

        if (tSys1.isTmp())
        {
            sys1.rename(name);
            return tSys1;
        }
        else
        {
            return tmp<linearSystem<STypeR,Type,MeshType>>
            (
                linearSystem<STypeR,Type,MeshType>::New
                (
                    name,
                    sys1.x()
                )
            );
        }
    }
};

template
<
    class STypeR,
    class SType12,
    class Type,
    class MeshType
>
class reuseSysTmpTmp<STypeR, STypeR, SType12, STypeR, Type, MeshType>
{
public:

    static tmp<linearSystem<STypeR,Type,MeshType>> New
    (
        const tmp<linearSystem<STypeR,Type,MeshType>>& tSys1,
        const tmp<linearSystem<STypeR,Type,MeshType>>& tSys2,
        const word name
    )
    {
        linearSystem<STypeR,Type,MeshType>& sys1 =
            const_cast<linearSystem<STypeR,Type,MeshType>&>(tSys1());

        linearSystem<STypeR,Type,MeshType>& sys2 =
            const_cast<linearSystem<STypeR,Type,MeshType>&>(tSys2());

        if (tSys1.isTmp())
        {
            sys1.rename(name);
            return tSys1;
        }
        else if (tSys2.isTmp())
        {
            sys2.rename(name);
            return tSys2;
        }
        else
        {
            return tmp<linearSystem<STypeR,Type,MeshType>>
            (
                linearSystem<STypeR,Type,MeshType>::New
                (
                    name,
                    sys1.x()
                )
            );
        }
    }
};

}

}

}

#endif
