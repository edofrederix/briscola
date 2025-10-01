#include "arguments.H"

#include "faceVector.H"

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

    faceVector h1;
    faceVector h2(Zero);
    faceVector h3(v1,v2,v3,v4,v5,v6);
    faceVector h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    faceVector h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    faceVector h6(is);

    rm(fileName);

    if (h6[0] != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6[1] != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6[2] != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6[3] != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6[4] != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6[5] != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);

    if (h6.left()   != v1) FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    if (h6.right()  != v2) FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    if (h6.bottom() != v3) FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    if (h6.top()    != v4) FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    if (h6.aft()    != v5) FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    if (h6.fore()   != v6) FatalErrorInFunction << "test 1f failed" << abort(FatalError);

    vector v7, v8, v9, v10, v11, v12;

    v7 = h5[0];
    v8 = h5[1];
    v9 = h5[2];
    v10 = h5[3];
    v11 = h5[4];
    v12 = h5[5];

    if (v7 != v1) FatalErrorInFunction << "test 2a failed" << abort(FatalError);
    if (v8 != v2) FatalErrorInFunction << "test 2b failed" << abort(FatalError);
    if (v9 != v3) FatalErrorInFunction << "test 2c failed" << abort(FatalError);
    if (v10 != v4) FatalErrorInFunction << "test 2d failed" << abort(FatalError);
    if (v11 != v5) FatalErrorInFunction << "test 2e failed" << abort(FatalError);
    if (v12 != v6) FatalErrorInFunction << "test 2f failed" << abort(FatalError);

    h5[0] = v6;
    h5[1] = v5;
    h5[2] = v4;
    h5[3] = v3;
    h5[4] = v2;
    h5[5] = v1;

    if (h5[0] != v6) FatalErrorInFunction << "test 3a failed" << abort(FatalError);
    if (h5[1] != v5) FatalErrorInFunction << "test 3b failed" << abort(FatalError);
    if (h5[2] != v4) FatalErrorInFunction << "test 3c failed" << abort(FatalError);
    if (h5[3] != v3) FatalErrorInFunction << "test 3d failed" << abort(FatalError);
    if (h5[4] != v2) FatalErrorInFunction << "test 3e failed" << abort(FatalError);
    if (h5[5] != v1) FatalErrorInFunction << "test 3f failed" << abort(FatalError);

    if (h6.left()    != h6[0]) FatalErrorInFunction << "test 5a failed" << abort(FatalError);
    if (h6.right()   != h6[1]) FatalErrorInFunction << "test 5b failed" << abort(FatalError);
    if (h6.bottom()  != h6[2]) FatalErrorInFunction << "test 5c failed" << abort(FatalError);
    if (h6.top()     != h6[3]) FatalErrorInFunction << "test 5d failed" << abort(FatalError);
    if (h6.aft()     != h6[4]) FatalErrorInFunction << "test 5e failed" << abort(FatalError);
    if (h6.fore()    != h6[5]) FatalErrorInFunction << "test 5f failed" << abort(FatalError);

    faceVector h7(faceVector::uniform(v1));

    if (h7.left()    != v1) FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    if (h7.right()   != v1) FatalErrorInFunction << "test 4b failed" << abort(FatalError);
    if (h7.bottom()  != v1) FatalErrorInFunction << "test 4c failed" << abort(FatalError);
    if (h7.top()     != v1) FatalErrorInFunction << "test 4d failed" << abort(FatalError);
    if (h7.aft()     != v1) FatalErrorInFunction << "test 4e failed" << abort(FatalError);
    if (h7.fore()    != v1) FatalErrorInFunction << "test 4f failed" << abort(FatalError);

    if (h6.lower().x() != v1) FatalErrorInFunction << "test 5a failed" << abort(FatalError);
    if (h6.lower().y() != v3) FatalErrorInFunction << "test 5b failed" << abort(FatalError);
    if (h6.lower().z() != v5) FatalErrorInFunction << "test 5c failed" << abort(FatalError);
    if (h6.upper().x() != v2) FatalErrorInFunction << "test 5d failed" << abort(FatalError);
    if (h6.upper().y() != v4) FatalErrorInFunction << "test 5e failed" << abort(FatalError);
    if (h6.upper().z() != v6) FatalErrorInFunction << "test 5f failed" << abort(FatalError);

    h7 += h6;
    h7 -= h6;

    h7 = Zero;
    h7 = pTraits<faceVector>::one;
    h7 = pTraits<typename faceVector::cmpt>::one;

    h7 *= scalar(2.0);
    h7 /= scalar(2.0);

    h3 + h4;
    h3 - h4;

    2.0*h3;
    h3/2.0;

    h3 & h4;
    h3 && h4;
    h3 ^ h4;

    sphericalTensor t1(2);
    diagTensor t2(1,2,3);
    symmTensor t3(1,2,3,4,5,6);
    tensor t4(1,2,3,4,5,6,7,8,9);

    v1 * h3;

    v1 & h3;
    t1 & h3;
    t2 & h3;
    t3 & h3;
    t4 & h3;

    v1 && h3;

    h3 & v1;
    h3 & t1;
    h3 & t2;
    h3 & t3;
    h3 & t4;

    h3 && v1;

    v1 + h3;
    h3 + v1;

    v1 - h3;
    h3 - v1;

    v1 & h3;
    h3 & v1;

    v1 ^ h3;
    h3 ^ v1;

    magSqr(h3);
    mag(h3);
    cmptMultiply(h3,h4);
    cmptPow(h3,h4);
    cmptDivide(h3,h4);
    cmptMax(h3);
    cmptMin(h3);
    cmptSum(h3);
    cmptAv(h3);
    cmptProduct(h3);
    cmptSqr(h3);
    cmptMag(h3);
    max(h3,h4);
    min(h3,h4);
    minMod(h3,h4);
}
