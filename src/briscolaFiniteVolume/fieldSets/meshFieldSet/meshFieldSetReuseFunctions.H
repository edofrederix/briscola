#ifndef meshFieldSetReuseFunctions_H
#define meshFieldSetReuseFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class reuseSetTmp
{
public:

    static tmp<Type> New
    (
        const tmp<Type>& tSet1,
        const word name
    )
    {
        Type& set1 = const_cast<Type&>(tSet1());

        if (tSet1.isTmp())
        {
            set1.rename(name);
            return tSet1;
        }
        else
        {
            return Type::New(name, set1.fvMsh());
        }
    }
};

template<class Type>
class reuseSetTmpTmp
{
public:

    static tmp<Type> New
    (
        const tmp<Type>& tSet1,
        const tmp<Type>& tSet2,
        const word name
    )
    {
        Type& set1 = const_cast<Type&>(tSet1());
        Type& set2 = const_cast<Type&>(tSet2());

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
            return Type::New(name, set1.fvMsh());
        }
    }
};

}

}

}

#endif
