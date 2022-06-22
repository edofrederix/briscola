#ifndef blockM_H
#define blockM_H

#include "error.H"

namespace Foam
{

namespace briscola
{

#ifdef FULLDEBUG

template<class Type1, class Type2>
void checkMatrices
(
    const block<Type1>& f1,
    const block<Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " block<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and block<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << endl << " for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3>
void checkMatrices
(
    const block<Type1>& f1,
    const block<Type2>& f2,
    const block<Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorInFunction
            << " block<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", block<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and block<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << endl << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<class Type1, class Type2>
void checkMatrices
(
    const block<Type1>&,
    const block<Type2>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3>
void checkMatrices
(
    const block<Type1>&,
    const block<Type2>&,
    const block<Type3>&,
    const char*
)
{}

#endif

}

}

#endif
