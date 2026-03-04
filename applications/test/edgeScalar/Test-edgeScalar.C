#include "arguments.H"

#include "edgeScalar.H"

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
    const scalar v9(9);
    const scalar v10(10);
    const scalar v11(11);
    const scalar v12(12);

    edgeScalar h1;
    edgeScalar h2(Zero);
    edgeScalar h3(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12);
    edgeScalar h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    edgeScalar h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    edgeScalar h6(is);

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

    pow3(h3);
    pow4(h3);
    pow5(h3);
    pow6(h3);
    pow025(h3);
    sqrt(h3);
    cbrt(h3);
    sign(h3);
    pos(h3);
    pos0(h3);
    neg(h3);
    neg0(h3);
    posPart(h3);
    negPart(h3);
    exp(h3);
    log(h3);
    log10(h3);
    sin(h3);
    cos(h3);
    tan(h3);
    sinh(h3);
    cosh(h3);
    tanh(h3);
    erf(h3);
    erfc(h3);
    lgamma(h3);
    j0(h3);
    j1(h3);
    y0(h3);
    y1(h3);

    if (!sigFpeEnabled())
    {
        asin(h3);
        acos(h3);
        atan(h3);
        asinh(h3);
        acosh(h3);
        atanh(h3);
    }
}
