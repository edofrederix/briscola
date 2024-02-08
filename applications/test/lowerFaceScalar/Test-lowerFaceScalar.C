#include "arguments.H"

#include "lowerFaceScalar.H"

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
#include "faceScalar.H"

using namespace Foam;
using namespace briscola;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    const scalar v1(1);
    const scalar v2(2);
    const scalar v3(3);

    lowerFaceScalar h1;
    lowerFaceScalar h2(Zero);
    lowerFaceScalar h3(v1,v2,v3);
    lowerFaceScalar h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    lowerFaceScalar h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    lowerFaceScalar h6(is);

    rm(fileName);

    if (h6.component(0) != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.component(1) != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.component(2) != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    if (h6.left()   != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.bottom() != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.aft()    != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    scalar v4, v5, v6;

    h5.component(v4, 0);
    h5.component(v5, 1);
    h5.component(v6, 2);

    if (v4 != v1) FatalErrorInFunction << "test 2a failed" << abort(FatalError);
    if (v5 != v2) FatalErrorInFunction << "test 2b failed" << abort(FatalError);
    if (v6 != v3) FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    h5.replace(0, v3);
    h5.replace(1, v2);
    h5.replace(2, v1);

    if (h5.component(0) != v3) FatalErrorInFunction << "test 3a failed" << abort(FatalError);
    if (h5.component(1) != v2) FatalErrorInFunction << "test 3b failed" << abort(FatalError);
    if (h5.component(2) != v1) FatalErrorInFunction << "test 3c failed" << abort(FatalError);

    if (h6.left()    != h6[0]) FatalErrorInFunction << "test 5a failed" << abort(FatalError);
    if (h6.bottom()  != h6[1]) FatalErrorInFunction << "test 5b failed" << abort(FatalError);
    if (h6.aft()     != h6[2]) FatalErrorInFunction << "test 5c failed" << abort(FatalError);

    lowerFaceScalar h7(lowerFaceScalar::uniform(v1));

    if (h7.left()    != v1) FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    if (h7.bottom()  != v1) FatalErrorInFunction << "test 4b failed" << abort(FatalError);
    if (h7.aft()     != v1) FatalErrorInFunction << "test 4c failed" << abort(FatalError);

    h7 += h6;

    if (h7.left()    != v1+v1) FatalErrorInFunction << "test 5d failed" << abort(FatalError);
    if (h7.bottom()  != v1+v2) FatalErrorInFunction << "test 5e failed" << abort(FatalError);
    if (h7.aft()     != v1+v3) FatalErrorInFunction << "test 5f failed" << abort(FatalError);

    h7 -= h6;

    if (h7.left()    != v1) FatalErrorInFunction << "test 6a failed" << abort(FatalError);
    if (h7.bottom()  != v1) FatalErrorInFunction << "test 6b failed" << abort(FatalError);
    if (h7.aft()     != v1) FatalErrorInFunction << "test 6c failed" << abort(FatalError);

    h7 = Zero;

    if (h7.left()    != 0) FatalErrorInFunction << "test 7a failed" << abort(FatalError);
    if (h7.bottom()  != 0) FatalErrorInFunction << "test 7b failed" << abort(FatalError);
    if (h7.aft()     != 0) FatalErrorInFunction << "test 7c failed" << abort(FatalError);

    h7 = pTraits<lowerFaceScalar>::one;

    if (h7.left()    != 1) FatalErrorInFunction << "test 8a failed" << abort(FatalError);
    if (h7.bottom()  != 1) FatalErrorInFunction << "test 8b failed" << abort(FatalError);
    if (h7.aft()     != 1) FatalErrorInFunction << "test 8c failed" << abort(FatalError);

    h7 = pTraits<typename lowerFaceScalar::cmpt>::one;

    if (h7.left()    != 1) FatalErrorInFunction << "test 8a failed" << abort(FatalError);
    if (h7.bottom()  != 1) FatalErrorInFunction << "test 8b failed" << abort(FatalError);
    if (h7.aft()     != 1) FatalErrorInFunction << "test 8c failed" << abort(FatalError);

    h7 = h6;
    h7 *= scalar(2.0);

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 9a failed" << abort(FatalError);
    if (h7.bottom()  != v2*2) FatalErrorInFunction << "test 9b failed" << abort(FatalError);
    if (h7.aft()     != v3*2) FatalErrorInFunction << "test 9c failed" << abort(FatalError);

    h7 /= scalar(2.0);

    if (h7.left()    != v1) FatalErrorInFunction << "test 10a failed" << abort(FatalError);
    if (h7.bottom()  != v2) FatalErrorInFunction << "test 10b failed" << abort(FatalError);
    if (h7.aft()     != v3) FatalErrorInFunction << "test 10c failed" << abort(FatalError);

    h7 = h6 + h7;

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 10a failed" << abort(FatalError);
    if (h7.bottom()  != v2*2) FatalErrorInFunction << "test 10b failed" << abort(FatalError);
    if (h7.aft()     != v3*2) FatalErrorInFunction << "test 10c failed" << abort(FatalError);

    h7 = h7 - h6;

    if (h7.left()    != v1) FatalErrorInFunction << "test 11a failed" << abort(FatalError);
    if (h7.bottom()  != v2) FatalErrorInFunction << "test 11b failed" << abort(FatalError);
    if (h7.aft()     != v3) FatalErrorInFunction << "test 11c failed" << abort(FatalError);

    h7 = 2.0*h7;

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 11a failed" << abort(FatalError);
    if (h7.bottom()  != v2*2) FatalErrorInFunction << "test 11b failed" << abort(FatalError);
    if (h7.aft()     != v3*2) FatalErrorInFunction << "test 11c failed" << abort(FatalError);

    h7 = h7/2.0;

    if (h7.left()    != v1) FatalErrorInFunction << "test 12a failed" << abort(FatalError);
    if (h7.bottom()  != v2) FatalErrorInFunction << "test 12b failed" << abort(FatalError);
    if (h7.aft()     != v3) FatalErrorInFunction << "test 12c failed" << abort(FatalError);

    faceScalar f1(v1,v2,v3,v4,v5,v6);

    h7 = f1;

    if (h7.left()    != v1) FatalErrorInFunction << "test 13a failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 13b failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 13c failed" << abort(FatalError);

    h7 += f1;

    if (h7.left()    != 2*v1) FatalErrorInFunction << "test 14a failed" << abort(FatalError);
    if (h7.bottom()  != 2*v3) FatalErrorInFunction << "test 14b failed" << abort(FatalError);
    if (h7.aft()     != 2*v5) FatalErrorInFunction << "test 14c failed" << abort(FatalError);

    h7 -= f1;

    if (h7.left()    != v1) FatalErrorInFunction << "test 15a failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 15b failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 15c failed" << abort(FatalError);

    h7 = h7*f1;

    if (h7.left()    != sqr(v1)) FatalErrorInFunction << "test 16a failed" << abort(FatalError);
    if (h7.bottom()  != sqr(v3)) FatalErrorInFunction << "test 16b failed" << abort(FatalError);
    if (h7.aft()     != sqr(v5)) FatalErrorInFunction << "test 16c failed" << abort(FatalError);

    h7 = h7/f1;

    if (h7.left()    != v1) FatalErrorInFunction << "test 17a failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 17b failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 17c failed" << abort(FatalError);

    h7 = cmptMultiply(h7,f1);

    if (h7.left()    != sqr(v1)) FatalErrorInFunction << "test 17a failed" << abort(FatalError);
    if (h7.bottom()  != sqr(v3)) FatalErrorInFunction << "test 17b failed" << abort(FatalError);
    if (h7.aft()     != sqr(v5)) FatalErrorInFunction << "test 17c failed" << abort(FatalError);

    h7 = cmptDivide(h7,f1);

    if (h7.left()    != v1) FatalErrorInFunction << "test 18a failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 18b failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 18c failed" << abort(FatalError);

    h7+f1;
    f1+h7;

    h7-f1;
    f1-h7;

    h7*f1;
    f1*h7;

    h7/f1;
    f1/h7;

    cmptPow(h7,f1);
    cmptPow(f1,h7);

    max(h7,f1);
    max(f1,h7);

    min(h7,f1);
    min(f1,h7);

    minMod(h7,f1);
    minMod(f1,h7);

    h7 == f1;
    f1 == h7;

    h7 != f1;
    f1 != h7;

    h7 > f1;
    f1 > h7;

    h7 < f1;
    f1 < h7;

    h7 >= f1;
    f1 >= h7;

    h7 <= f1;
    f1 <= h7;
}
