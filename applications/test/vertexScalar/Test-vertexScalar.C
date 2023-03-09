#include "arguments.H"

#include "vertexScalar.H"

#include "fileOperation.H"
#include "OSspecific.H"

#include "IFstream.H"
#include "OFstream.H"

#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "diagTensor.H"

using namespace Foam;
using namespace briscola;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    const scalar v1(1);
    const scalar v2(2);
    const scalar v3(3);
    const scalar v4(4);
    const scalar v5(5);
    const scalar v6(6);
    const scalar v7(7);
    const scalar v8(8);

    vertexScalar h1;
    vertexScalar h2(Zero);
    vertexScalar h3(v1,v2,v3,v4,v5,v6,v7,v8);
    vertexScalar h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    vertexScalar h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    vertexScalar h6(is);

    rm(fileName);

    if (h6.component(0) != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.component(1) != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.component(2) != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.component(3) != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.component(4) != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.component(5) != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6.component(6) != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6.component(7) != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);

    if (h6.v0() != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.v1() != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.v2() != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.v3() != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.v4() != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.v5() != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6.v6() != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6.v7() != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);
}
