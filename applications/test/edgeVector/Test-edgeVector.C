#include "arguments.H"

#include "edgeVector.H"

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
    const vector v9(3,4,5);
    const vector v10(4,5,6);
    const vector v11(5,6,7);
    const vector v12(6,7,8);

    edgeVector h1;
    edgeVector h2(Zero);
    edgeVector h3(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12);
    edgeVector h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    edgeVector h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    edgeVector h6(is);

    rm(fileName);

    if (h6[0] != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6[1] != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6[2] != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6[3] != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6[4] != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6[5] != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6[6] != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6[7] != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);
    if (h6[8] != v9) FatalErrorInFunction << "test 1i failed" << abort(FatalError);
    if (h6[9] != v10) FatalErrorInFunction << "test 1j failed" << abort(FatalError);
    if (h6[10] != v11) FatalErrorInFunction << "test 1k failed" << abort(FatalError);
    if (h6[11] != v12) FatalErrorInFunction << "test 1l failed" << abort(FatalError);

    if (h6.ba() != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.ta() != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.bf() != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.tf() != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.la() != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.ra() != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    if (h6.lf() != v7) FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    if (h6.rf() != v8) FatalErrorInFunction << "test 1h failed" << abort(FatalError);
    if (h6.lb() != v9) FatalErrorInFunction << "test 1i failed" << abort(FatalError);
    if (h6.rb() != v10) FatalErrorInFunction << "test 1j failed" << abort(FatalError);
    if (h6.lt() != v11) FatalErrorInFunction << "test 1k failed" << abort(FatalError);
    if (h6.rt() != v12) FatalErrorInFunction << "test 1l failed" << abort(FatalError);
}
