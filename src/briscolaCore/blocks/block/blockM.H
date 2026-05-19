#ifndef blockM_H
#define blockM_H

#include "error.H"

namespace Foam
{

namespace briscola
{

#ifdef FULLDEBUG

template<class Type1, class Type2, int P>
void checkBlocks
(
    const block<Type1,P>& f1,
    const block<Type2,P>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
        FatalErrorInFunction
            << " block<"<<pTraits<Type1>::typeName
            << "," << Foam::name(P) << "> f1("<<f1.size()<<')'
            << " and block<"<<pTraits<Type2>::typeName
            << "," << Foam::name(P) << "> f2("<<f2.size()<<')'
            << endl << " for operation " << op
            << abort(FatalError);
}

template<class Type1, class Type2, class Type3, int P>
void checkBlocks
(
    const block<Type1,P>& f1,
    const block<Type2,P>& f2,
    const block<Type3,P>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
        FatalErrorInFunction
            << " block<"<<pTraits<Type1>::typeName
            < <"," << Foam::name(P) << "> f1("<<f1.size()<<')'
            << ", block<"<<pTraits<Type2>::typeName
            < <"," << Foam::name(P) << "> f2("<<f2.size()<<')'
            << " and block<"<<pTraits<Type3>::typeName
            < <"," << Foam::name(P) << "> f3("<<f3.size()<<')'
            << endl << "    for operation " << op
            << abort(FatalError);
}

#else

template<class Type1, class Type2, int P>
void checkBlocks
(
    const block<Type1,P>&,
    const block<Type2,P>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3, int P>
void checkBlocks
(
    const block<Type1,P>&,
    const block<Type2,P>&,
    const block<Type3,P>&,
    const char*
)
{}

#endif

}

}

#endif
