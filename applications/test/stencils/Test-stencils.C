#include "arguments.H"

#include "stencil.H"
#include "symmStencil.H"
#include "diagStencil.H"

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

    st2 += st1;
    st2 -= st1;

    st2 = Zero;
    st2 = st1;
    st2 = 2.0;

    st2 *= 2.0;
    st2 /= 2.0;
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    testConstructors<stencil>();
    testConstructors<symmStencil>();
    testConstructors<diagStencil>();

    stencil st1(1,2,3,4,5,6,7);
    symmStencil st2(1,2,3,4);
    diagStencil st3(1);

    testMemberOperators<stencil>(st1);
    testMemberOperators<symmStencil>(st2);
    testMemberOperators<diagStencil>(st3);

    // Stencil tests

    if (st1.center() != 1) FatalErrorInFunction << "test 2a failed" << abort(FatalError);
    if (st1.left()   != 2) FatalErrorInFunction << "test 2b failed" << abort(FatalError);
    if (st1.right()  != 3) FatalErrorInFunction << "test 2c failed" << abort(FatalError);
    if (st1.bottom() != 4) FatalErrorInFunction << "test 2d failed" << abort(FatalError);
    if (st1.top()    != 5) FatalErrorInFunction << "test 2e failed" << abort(FatalError);
    if (st1.aft()    != 6) FatalErrorInFunction << "test 2f failed" << abort(FatalError);
    if (st1.fore()   != 7) FatalErrorInFunction << "test 2g failed" << abort(FatalError);

    if (st1.center() != st1[0]) FatalErrorInFunction << "test 3a failed" << abort(FatalError);
    if (st1.left()   != st1[1]) FatalErrorInFunction << "test 3b failed" << abort(FatalError);
    if (st1.right()  != st1[2]) FatalErrorInFunction << "test 3c failed" << abort(FatalError);
    if (st1.bottom() != st1[3]) FatalErrorInFunction << "test 3d failed" << abort(FatalError);
    if (st1.top()    != st1[4]) FatalErrorInFunction << "test 3e failed" << abort(FatalError);
    if (st1.aft()    != st1[5]) FatalErrorInFunction << "test 3f failed" << abort(FatalError);
    if (st1.fore()   != st1[6]) FatalErrorInFunction << "test 3g failed" << abort(FatalError);

    st1 += st2;
    st1 += st3;

    st1 -= st2;
    st1 -= st3;

    st1 = st2;
    st1 = st3;

    // SymmStencil tests

    if (st2.center() != 1) FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    if (st2.left()   != 2) FatalErrorInFunction << "test 4b failed" << abort(FatalError);
    if (st2.bottom() != 3) FatalErrorInFunction << "test 4c failed" << abort(FatalError);
    if (st2.aft()    != 4) FatalErrorInFunction << "test 4d failed" << abort(FatalError);

    if (st2.center() != st2[0]) FatalErrorInFunction << "test 5a failed" << abort(FatalError);
    if (st2.left()   != st2[1]) FatalErrorInFunction << "test 5b failed" << abort(FatalError);
    if (st2.bottom() != st2[2]) FatalErrorInFunction << "test 5c failed" << abort(FatalError);
    if (st2.aft()    != st2[3]) FatalErrorInFunction << "test 5d failed" << abort(FatalError);

    st2 += st3;
    st2 -= st3;
    st2 = st3;

    // DiagStencil tests

    if (st3.center() != 1) FatalErrorInFunction << "test 6 failed" << abort(FatalError);
    if (st3.center() != st3[0]) FatalErrorInFunction << "test 7 failed" << abort(FatalError);
}
