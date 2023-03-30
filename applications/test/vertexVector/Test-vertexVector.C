#include "arguments.H"

#include "vertexVector.H"

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

    const vector v1(1,2,3);
    const vector v2(2,3,4);
    const vector v3(3,4,5);
    const vector v4(4,5,6);
    const vector v5(5,6,7);
    const vector v6(6,7,8);
    const vector v7(1,2,3);
    const vector v8(2,3,4);

    vertexVector h1;
    vertexVector h2(Zero);
    vertexVector h3(v1,v2,v3,v4,v5,v6,v7,v8);
    vertexVector h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    vertexVector h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    vertexVector h6(is);

    rm(fileName);

    if (h6.component(0) != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.component(1) != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.component(2) != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.component(3) != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.component(4) != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.component(5) != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6.component(6) != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6.component(7) != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);

    if (h6.lba() != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.rba() != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.lta() != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.rta() != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.lbf() != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.rbf() != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6.ltf() != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6.rtf() != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);
}
