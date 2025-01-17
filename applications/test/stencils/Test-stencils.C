#include "arguments.H"

#include "stencil.H"
#include "diagStencil.H"
#include "faceScalar.H"

#include "fileOperation.H"
#include "OSspecific.H"

#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;
using namespace briscola;

template<class Type>
void testConstructors()
{
    Type st1;
    Type st2(Zero);
    Type st3(2.0);

    const word fileName =
        "dummy-"
      + word(pTraits<Type>::typeName)
      + "-"
      + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << st3 << endl;

    IFstream is(fileName);

    Type st4(is);

    if (st3 != st4) FatalErrorInFunction << "test 1 failed" << abort(FatalError);

    rm(fileName);

    Type st5(Type::uniform(3.0));
}

template<class Type>
void testMemberOperators(const Type& st1)
{
    Type st2(2.0*st1);

    st2 = pTraits<Type>::one;
    st2 = pTraits<typename Type::cmpt>::one;

    st2 += st1;
    st2 -= st1;

    st2 = Zero;
    st2 = st1;
    st2 = 2.0;

    st2 += 2.0;
    st2 *= 2.0;
    st2 /= 2.0;

    st2 + 2.0;
    st2 - 2.0;

    st1 * 2.0;
    st2 / 2.0;
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    testConstructors<stencil>();
    testConstructors<diagStencil>();

    stencil st2(1,2,3,4,5,6,7);
    diagStencil st3(1);

    testMemberOperators<stencil>(st2);
    testMemberOperators<diagStencil>(st3);

    // Stencil tests

    if (st2.center() != 1) FatalErrorInFunction << "test 2a failed" << abort(FatalError);
    if (st2.left()   != 2) FatalErrorInFunction << "test 2b failed" << abort(FatalError);
    if (st2.right()  != 3) FatalErrorInFunction << "test 2c failed" << abort(FatalError);
    if (st2.bottom() != 4) FatalErrorInFunction << "test 2d failed" << abort(FatalError);
    if (st2.top()    != 5) FatalErrorInFunction << "test 2e failed" << abort(FatalError);
    if (st2.aft()    != 6) FatalErrorInFunction << "test 2f failed" << abort(FatalError);
    if (st2.fore()   != 7) FatalErrorInFunction << "test 2g failed" << abort(FatalError);

    if (st2.center() != st2[0]) FatalErrorInFunction << "test 3a failed" << abort(FatalError);
    if (st2.left()   != st2[1]) FatalErrorInFunction << "test 3b failed" << abort(FatalError);
    if (st2.right()  != st2[2]) FatalErrorInFunction << "test 3c failed" << abort(FatalError);
    if (st2.bottom() != st2[3]) FatalErrorInFunction << "test 3d failed" << abort(FatalError);
    if (st2.top()    != st2[4]) FatalErrorInFunction << "test 3e failed" << abort(FatalError);
    if (st2.aft()    != st2[5]) FatalErrorInFunction << "test 3f failed" << abort(FatalError);
    if (st2.fore()   != st2[6]) FatalErrorInFunction << "test 3g failed" << abort(FatalError);

    st2 += st3;
    st2 -= st3;
    st2 = st3;

    faceScalar hs(1,2,3,4,5,6);

    st2 += hs;
    st2 -= hs;
    st2 = hs;

    st2 + hs;
    st2 - hs;
    st2 * hs;
}
