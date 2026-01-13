#include "arguments.H"

#include "faceScalar.H"

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

    faceScalar h1;
    faceScalar h2(Zero);
    faceScalar h3(v1,v2,v3,v4,v5,v6);
    faceScalar h4(h3);

    const word fileName = "dummy-" + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << h4 << endl;

    faceScalar h5;

    IFstream is(fileName);

    is >> h5;

    is.rewind();

    faceScalar h6(is);

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

    scalar v7, v8, v9, v10, v11, v12;

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

    faceScalar h7(faceScalar::uniform(v1));

    if (h7.left()    != v1) FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    if (h7.right()   != v1) FatalErrorInFunction << "test 4b failed" << abort(FatalError);
    if (h7.bottom()  != v1) FatalErrorInFunction << "test 4c failed" << abort(FatalError);
    if (h7.top()     != v1) FatalErrorInFunction << "test 4d failed" << abort(FatalError);
    if (h7.aft()     != v1) FatalErrorInFunction << "test 4e failed" << abort(FatalError);
    if (h7.fore()    != v1) FatalErrorInFunction << "test 4f failed" << abort(FatalError);

    h7 += h6;

    if (h7.left()    != v1+v1) FatalErrorInFunction << "test 5a failed" << abort(FatalError);
    if (h7.right()   != v1+v2) FatalErrorInFunction << "test 5b failed" << abort(FatalError);
    if (h7.bottom()  != v1+v3) FatalErrorInFunction << "test 5c failed" << abort(FatalError);
    if (h7.top()     != v1+v4) FatalErrorInFunction << "test 5d failed" << abort(FatalError);
    if (h7.aft()     != v1+v5) FatalErrorInFunction << "test 5e failed" << abort(FatalError);
    if (h7.fore()    != v1+v6) FatalErrorInFunction << "test 5f failed" << abort(FatalError);

    h7 -= h6;

    if (h7.left()    != v1) FatalErrorInFunction << "test 6a failed" << abort(FatalError);
    if (h7.right()   != v1) FatalErrorInFunction << "test 6b failed" << abort(FatalError);
    if (h7.bottom()  != v1) FatalErrorInFunction << "test 6c failed" << abort(FatalError);
    if (h7.top()     != v1) FatalErrorInFunction << "test 6d failed" << abort(FatalError);
    if (h7.aft()     != v1) FatalErrorInFunction << "test 6e failed" << abort(FatalError);
    if (h7.fore()    != v1) FatalErrorInFunction << "test 6f failed" << abort(FatalError);

    h7 = Zero;

    if (h7.left()    != 0) FatalErrorInFunction << "test 7a failed" << abort(FatalError);
    if (h7.right()   != 0) FatalErrorInFunction << "test 7b failed" << abort(FatalError);
    if (h7.bottom()  != 0) FatalErrorInFunction << "test 7c failed" << abort(FatalError);
    if (h7.top()     != 0) FatalErrorInFunction << "test 7d failed" << abort(FatalError);
    if (h7.aft()     != 0) FatalErrorInFunction << "test 7e failed" << abort(FatalError);
    if (h7.fore()    != 0) FatalErrorInFunction << "test 7f failed" << abort(FatalError);

    h7 = pTraits<faceScalar>::one;

    if (h7.left()    != 1) FatalErrorInFunction << "test 8a failed" << abort(FatalError);
    if (h7.right()   != 1) FatalErrorInFunction << "test 8b failed" << abort(FatalError);
    if (h7.bottom()  != 1) FatalErrorInFunction << "test 8c failed" << abort(FatalError);
    if (h7.top()     != 1) FatalErrorInFunction << "test 8d failed" << abort(FatalError);
    if (h7.aft()     != 1) FatalErrorInFunction << "test 8e failed" << abort(FatalError);
    if (h7.fore()    != 1) FatalErrorInFunction << "test 8f failed" << abort(FatalError);

    h7 = pTraits<typename faceScalar::cmpt>::one;

    if (h7.left()    != 1) FatalErrorInFunction << "test 8a failed" << abort(FatalError);
    if (h7.right()   != 1) FatalErrorInFunction << "test 8b failed" << abort(FatalError);
    if (h7.bottom()  != 1) FatalErrorInFunction << "test 8c failed" << abort(FatalError);
    if (h7.top()     != 1) FatalErrorInFunction << "test 8d failed" << abort(FatalError);
    if (h7.aft()     != 1) FatalErrorInFunction << "test 8e failed" << abort(FatalError);
    if (h7.fore()    != 1) FatalErrorInFunction << "test 8f failed" << abort(FatalError);

    h7 = h6;
    h7 *= scalar(2.0);

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 9a failed" << abort(FatalError);
    if (h7.right()   != v2*2) FatalErrorInFunction << "test 9b failed" << abort(FatalError);
    if (h7.bottom()  != v3*2) FatalErrorInFunction << "test 9c failed" << abort(FatalError);
    if (h7.top()     != v4*2) FatalErrorInFunction << "test 9d failed" << abort(FatalError);
    if (h7.aft()     != v5*2) FatalErrorInFunction << "test 9e failed" << abort(FatalError);
    if (h7.fore()    != v6*2) FatalErrorInFunction << "test 9f failed" << abort(FatalError);

    h7 /= scalar(2.0);

    if (h7.left()    != v1) FatalErrorInFunction << "test 10a failed" << abort(FatalError);
    if (h7.right()   != v2) FatalErrorInFunction << "test 10b failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 10c failed" << abort(FatalError);
    if (h7.top()     != v4) FatalErrorInFunction << "test 10d failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 10e failed" << abort(FatalError);
    if (h7.fore()    != v6) FatalErrorInFunction << "test 10f failed" << abort(FatalError);

    h7 = h6 + h7;

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 10a failed" << abort(FatalError);
    if (h7.right()   != v2*2) FatalErrorInFunction << "test 10b failed" << abort(FatalError);
    if (h7.bottom()  != v3*2) FatalErrorInFunction << "test 10c failed" << abort(FatalError);
    if (h7.top()     != v4*2) FatalErrorInFunction << "test 10d failed" << abort(FatalError);
    if (h7.aft()     != v5*2) FatalErrorInFunction << "test 10e failed" << abort(FatalError);
    if (h7.fore()    != v6*2) FatalErrorInFunction << "test 10f failed" << abort(FatalError);

    h7 = h7 - h6;

    if (h7.left()    != v1) FatalErrorInFunction << "test 11a failed" << abort(FatalError);
    if (h7.right()   != v2) FatalErrorInFunction << "test 11b failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 11c failed" << abort(FatalError);
    if (h7.top()     != v4) FatalErrorInFunction << "test 11d failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 11e failed" << abort(FatalError);
    if (h7.fore()    != v6) FatalErrorInFunction << "test 11f failed" << abort(FatalError);

    h7 = 2.0*h7;

    if (h7.left()    != v1*2) FatalErrorInFunction << "test 11a failed" << abort(FatalError);
    if (h7.right()   != v2*2) FatalErrorInFunction << "test 11b failed" << abort(FatalError);
    if (h7.bottom()  != v3*2) FatalErrorInFunction << "test 11c failed" << abort(FatalError);
    if (h7.top()     != v4*2) FatalErrorInFunction << "test 11d failed" << abort(FatalError);
    if (h7.aft()     != v5*2) FatalErrorInFunction << "test 11e failed" << abort(FatalError);
    if (h7.fore()    != v6*2) FatalErrorInFunction << "test 11f failed" << abort(FatalError);

    h7 = h7/2.0;

    if (h7.left()    != v1) FatalErrorInFunction << "test 12a failed" << abort(FatalError);
    if (h7.right()   != v2) FatalErrorInFunction << "test 12b failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 12c failed" << abort(FatalError);
    if (h7.top()     != v4) FatalErrorInFunction << "test 12d failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 12e failed" << abort(FatalError);
    if (h7.fore()    != v6) FatalErrorInFunction << "test 12f failed" << abort(FatalError);

    vector vv1(1,2,3);
    sphericalTensor t1(2);
    diagTensor t2(1,2,3);
    symmTensor t3(1,2,3,4,5,6);
    tensor t4(1,2,3,4,5,6,7,8,9);

    FaceSpace<vector> fv7 = vv1 * h7;

    if (fv7.left()    != v1*vv1) FatalErrorInFunction << "test 13a failed" << abort(FatalError);
    if (fv7.right()   != v2*vv1) FatalErrorInFunction << "test 13b failed" << abort(FatalError);
    if (fv7.bottom()  != v3*vv1) FatalErrorInFunction << "test 13c failed" << abort(FatalError);
    if (fv7.top()     != v4*vv1) FatalErrorInFunction << "test 13d failed" << abort(FatalError);
    if (fv7.aft()     != v5*vv1) FatalErrorInFunction << "test 13e failed" << abort(FatalError);
    if (fv7.fore()    != v6*vv1) FatalErrorInFunction << "test 13f failed" << abort(FatalError);

    FaceSpace<sphericalTensor> fst7 = t1 * h7;

    if (fst7.left()    != v1*t1) FatalErrorInFunction << "test 14a failed" << abort(FatalError);
    if (fst7.right()   != v2*t1) FatalErrorInFunction << "test 14b failed" << abort(FatalError);
    if (fst7.bottom()  != v3*t1) FatalErrorInFunction << "test 14c failed" << abort(FatalError);
    if (fst7.top()     != v4*t1) FatalErrorInFunction << "test 14d failed" << abort(FatalError);
    if (fst7.aft()     != v5*t1) FatalErrorInFunction << "test 14e failed" << abort(FatalError);
    if (fst7.fore()    != v6*t1) FatalErrorInFunction << "test 14f failed" << abort(FatalError);

    // Bug in OpenFOAM:
    // FaceSpace<diagTensor> fdt7 = t2 * h7;

    // if (fdt7.left()    != v1*t2) FatalErrorInFunction << "test 15a failed" << abort(FatalError);
    // if (fdt7.right()   != v2*t2) FatalErrorInFunction << "test 15b failed" << abort(FatalError);
    // if (fdt7.bottom()  != v3*t2) FatalErrorInFunction << "test 15c failed" << abort(FatalError);
    // if (fdt7.top()     != v4*t2) FatalErrorInFunction << "test 15d failed" << abort(FatalError);
    // if (fdt7.aft()     != v5*t2) FatalErrorInFunction << "test 15e failed" << abort(FatalError);
    // if (fdt7.fore()    != v6*t2) FatalErrorInFunction << "test 15f failed" << abort(FatalError);

    FaceSpace<symmTensor> fsmt7 = t3 * h7;

    if (fsmt7.left()    != v1*t3) FatalErrorInFunction << "test 16a failed" << abort(FatalError);
    if (fsmt7.right()   != v2*t3) FatalErrorInFunction << "test 16b failed" << abort(FatalError);
    if (fsmt7.bottom()  != v3*t3) FatalErrorInFunction << "test 16c failed" << abort(FatalError);
    if (fsmt7.top()     != v4*t3) FatalErrorInFunction << "test 16d failed" << abort(FatalError);
    if (fsmt7.aft()     != v5*t3) FatalErrorInFunction << "test 16e failed" << abort(FatalError);
    if (fsmt7.fore()    != v6*t3) FatalErrorInFunction << "test 16f failed" << abort(FatalError);

    FaceSpace<tensor> ft7 = t4 * h7;

    if (ft7.left()    != v1*t4) FatalErrorInFunction << "test 17a failed" << abort(FatalError);
    if (ft7.right()   != v2*t4) FatalErrorInFunction << "test 17b failed" << abort(FatalError);
    if (ft7.bottom()  != v3*t4) FatalErrorInFunction << "test 17c failed" << abort(FatalError);
    if (ft7.top()     != v4*t4) FatalErrorInFunction << "test 17d failed" << abort(FatalError);
    if (ft7.aft()     != v5*t4) FatalErrorInFunction << "test 17e failed" << abort(FatalError);
    if (ft7.fore()    != v6*t4) FatalErrorInFunction << "test 17f failed" << abort(FatalError);

    FaceSpace<vector> fv7i = vv1 / h7;

    if (fv7i.left()    != vv1/v1) FatalErrorInFunction << "test 18a failed" << abort(FatalError);
    if (fv7i.right()   != vv1/v2) FatalErrorInFunction << "test 18b failed" << abort(FatalError);
    if (fv7i.bottom()  != vv1/v3) FatalErrorInFunction << "test 18c failed" << abort(FatalError);
    if (fv7i.top()     != vv1/v4) FatalErrorInFunction << "test 18d failed" << abort(FatalError);
    if (fv7i.aft()     != vv1/v5) FatalErrorInFunction << "test 18e failed" << abort(FatalError);
    if (fv7i.fore()    != vv1/v6) FatalErrorInFunction << "test 18f failed" << abort(FatalError);

    FaceSpace<sphericalTensor> fst7i = t1 / h7;

    if (fst7i.left()    != t1/v1) FatalErrorInFunction << "test 19a failed" << abort(FatalError);
    if (fst7i.right()   != t1/v2) FatalErrorInFunction << "test 19b failed" << abort(FatalError);
    if (fst7i.bottom()  != t1/v3) FatalErrorInFunction << "test 19c failed" << abort(FatalError);
    if (fst7i.top()     != t1/v4) FatalErrorInFunction << "test 19d failed" << abort(FatalError);
    if (fst7i.aft()     != t1/v5) FatalErrorInFunction << "test 19e failed" << abort(FatalError);
    if (fst7i.fore()    != t1/v6) FatalErrorInFunction << "test 19f failed" << abort(FatalError);

    FaceSpace<diagTensor> fdt7i = t2 / h7;

    if (fdt7i.left()    != t2/v1) FatalErrorInFunction << "test 20a failed" << abort(FatalError);
    if (fdt7i.right()   != t2/v2) FatalErrorInFunction << "test 20b failed" << abort(FatalError);
    if (fdt7i.bottom()  != t2/v3) FatalErrorInFunction << "test 20c failed" << abort(FatalError);
    if (fdt7i.top()     != t2/v4) FatalErrorInFunction << "test 20d failed" << abort(FatalError);
    if (fdt7i.aft()     != t2/v5) FatalErrorInFunction << "test 20e failed" << abort(FatalError);
    if (fdt7i.fore()    != t2/v6) FatalErrorInFunction << "test 20f failed" << abort(FatalError);

    FaceSpace<symmTensor> fsmt7i = t3 / h7;

    if (fsmt7i.left()    != t3/v1) FatalErrorInFunction << "test 21a failed" << abort(FatalError);
    if (fsmt7i.right()   != t3/v2) FatalErrorInFunction << "test 21b failed" << abort(FatalError);
    if (fsmt7i.bottom()  != t3/v3) FatalErrorInFunction << "test 21c failed" << abort(FatalError);
    if (fsmt7i.top()     != t3/v4) FatalErrorInFunction << "test 21d failed" << abort(FatalError);
    if (fsmt7i.aft()     != t3/v5) FatalErrorInFunction << "test 21e failed" << abort(FatalError);
    if (fsmt7i.fore()    != t3/v6) FatalErrorInFunction << "test 21f failed" << abort(FatalError);

    FaceSpace<tensor> ft7i = t4 / h7;

    if (ft7i.left()    != t4/v1) FatalErrorInFunction << "test 22a failed" << abort(FatalError);
    if (ft7i.right()   != t4/v2) FatalErrorInFunction << "test 22b failed" << abort(FatalError);
    if (ft7i.bottom()  != t4/v3) FatalErrorInFunction << "test 22c failed" << abort(FatalError);
    if (ft7i.top()     != t4/v4) FatalErrorInFunction << "test 22d failed" << abort(FatalError);
    if (ft7i.aft()     != t4/v5) FatalErrorInFunction << "test 22e failed" << abort(FatalError);
    if (ft7i.fore()    != t4/v6) FatalErrorInFunction << "test 22f failed" << abort(FatalError);

    fv7 = h7 * vv1;

    if (fv7.left()    != v1*vv1) FatalErrorInFunction << "test 23a failed" << abort(FatalError);
    if (fv7.right()   != v2*vv1) FatalErrorInFunction << "test 23b failed" << abort(FatalError);
    if (fv7.bottom()  != v3*vv1) FatalErrorInFunction << "test 23c failed" << abort(FatalError);
    if (fv7.top()     != v4*vv1) FatalErrorInFunction << "test 23d failed" << abort(FatalError);
    if (fv7.aft()     != v5*vv1) FatalErrorInFunction << "test 23e failed" << abort(FatalError);
    if (fv7.fore()    != v6*vv1) FatalErrorInFunction << "test 23f failed" << abort(FatalError);

    fst7 = h7 * t1;

    if (fst7.left()    != v1*t1) FatalErrorInFunction << "test 24a failed" << abort(FatalError);
    if (fst7.right()   != v2*t1) FatalErrorInFunction << "test 24b failed" << abort(FatalError);
    if (fst7.bottom()  != v3*t1) FatalErrorInFunction << "test 24c failed" << abort(FatalError);
    if (fst7.top()     != v4*t1) FatalErrorInFunction << "test 24d failed" << abort(FatalError);
    if (fst7.aft()     != v5*t1) FatalErrorInFunction << "test 24e failed" << abort(FatalError);
    if (fst7.fore()    != v6*t1) FatalErrorInFunction << "test 24f failed" << abort(FatalError);

    // Bug in OpenFOAM:
    // fdt7 = h7 * t2;

    // if (fdt7.left()    != v1*t2) FatalErrorInFunction << "test 25a failed" << abort(FatalError);
    // if (fdt7.right()   != v2*t2) FatalErrorInFunction << "test 25b failed" << abort(FatalError);
    // if (fdt7.bottom()  != v3*t2) FatalErrorInFunction << "test 25c failed" << abort(FatalError);
    // if (fdt7.top()     != v4*t2) FatalErrorInFunction << "test 25d failed" << abort(FatalError);
    // if (fdt7.aft()     != v5*t2) FatalErrorInFunction << "test 25e failed" << abort(FatalError);
    // if (fdt7.fore()    != v6*t2) FatalErrorInFunction << "test 25f failed" << abort(FatalError);

    fsmt7 = h7 * t3;

    if (fsmt7.left()    != v1*t3) FatalErrorInFunction << "test 26a failed" << abort(FatalError);
    if (fsmt7.right()   != v2*t3) FatalErrorInFunction << "test 26b failed" << abort(FatalError);
    if (fsmt7.bottom()  != v3*t3) FatalErrorInFunction << "test 26c failed" << abort(FatalError);
    if (fsmt7.top()     != v4*t3) FatalErrorInFunction << "test 26d failed" << abort(FatalError);
    if (fsmt7.aft()     != v5*t3) FatalErrorInFunction << "test 26e failed" << abort(FatalError);
    if (fsmt7.fore()    != v6*t3) FatalErrorInFunction << "test 26f failed" << abort(FatalError);

    ft7 = h7 * t4;

    if (ft7.left()    != v1*t4) FatalErrorInFunction << "test 27a failed" << abort(FatalError);
    if (ft7.right()   != v2*t4) FatalErrorInFunction << "test 27b failed" << abort(FatalError);
    if (ft7.bottom()  != v3*t4) FatalErrorInFunction << "test 27c failed" << abort(FatalError);
    if (ft7.top()     != v4*t4) FatalErrorInFunction << "test 27d failed" << abort(FatalError);
    if (ft7.aft()     != v5*t4) FatalErrorInFunction << "test 27e failed" << abort(FatalError);
    if (ft7.fore()    != v6*t4) FatalErrorInFunction << "test 27f failed" << abort(FatalError);

    h7 = v1 + h7;

    if (h7.left()    != v1+v1) FatalErrorInFunction << "test 28a failed" << abort(FatalError);
    if (h7.right()   != v2+v1) FatalErrorInFunction << "test 28b failed" << abort(FatalError);
    if (h7.bottom()  != v3+v1) FatalErrorInFunction << "test 28c failed" << abort(FatalError);
    if (h7.top()     != v4+v1) FatalErrorInFunction << "test 28d failed" << abort(FatalError);
    if (h7.aft()     != v5+v1) FatalErrorInFunction << "test 28e failed" << abort(FatalError);
    if (h7.fore()    != v6+v1) FatalErrorInFunction << "test 28f failed" << abort(FatalError);

    h7 = h7 + v1;

    if (h7.left()    != v1+v1*2) FatalErrorInFunction << "test 29a failed" << abort(FatalError);
    if (h7.right()   != v2+v1*2) FatalErrorInFunction << "test 29b failed" << abort(FatalError);
    if (h7.bottom()  != v3+v1*2) FatalErrorInFunction << "test 29c failed" << abort(FatalError);
    if (h7.top()     != v4+v1*2) FatalErrorInFunction << "test 29d failed" << abort(FatalError);
    if (h7.aft()     != v5+v1*2) FatalErrorInFunction << "test 29e failed" << abort(FatalError);
    if (h7.fore()    != v6+v1*2) FatalErrorInFunction << "test 29f failed" << abort(FatalError);

    h7 = v1 - h7;

    if (h7.left()    != v1-(v1+v1*2)) FatalErrorInFunction << "test 30a failed" << abort(FatalError);
    if (h7.right()   != v1-(v2+v1*2)) FatalErrorInFunction << "test 30b failed" << abort(FatalError);
    if (h7.bottom()  != v1-(v3+v1*2)) FatalErrorInFunction << "test 30c failed" << abort(FatalError);
    if (h7.top()     != v1-(v4+v1*2)) FatalErrorInFunction << "test 30d failed" << abort(FatalError);
    if (h7.aft()     != v1-(v5+v1*2)) FatalErrorInFunction << "test 30e failed" << abort(FatalError);
    if (h7.fore()    != v1-(v6+v1*2)) FatalErrorInFunction << "test 30f failed" << abort(FatalError);

    h7 = h7 - v1;

    if (h7.left()    != v1-(v1+v1*2)-v1) FatalErrorInFunction << "test 31a failed" << abort(FatalError);
    if (h7.right()   != v1-(v2+v1*2)-v1) FatalErrorInFunction << "test 31b failed" << abort(FatalError);
    if (h7.bottom()  != v1-(v3+v1*2)-v1) FatalErrorInFunction << "test 31c failed" << abort(FatalError);
    if (h7.top()     != v1-(v4+v1*2)-v1) FatalErrorInFunction << "test 31d failed" << abort(FatalError);
    if (h7.aft()     != v1-(v5+v1*2)-v1) FatalErrorInFunction << "test 31e failed" << abort(FatalError);
    if (h7.fore()    != v1-(v6+v1*2)-v1) FatalErrorInFunction << "test 31f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = h7 && h7;

    if (h7.left()    != (v1 && v1)) FatalErrorInFunction << "test 32a failed" << abort(FatalError);
    if (h7.right()   != (v2 && v2)) FatalErrorInFunction << "test 32b failed" << abort(FatalError);
    if (h7.bottom()  != (v3 && v3)) FatalErrorInFunction << "test 32c failed" << abort(FatalError);
    if (h7.top()     != (v4 && v4)) FatalErrorInFunction << "test 32d failed" << abort(FatalError);
    if (h7.aft()     != (v5 && v5)) FatalErrorInFunction << "test 32e failed" << abort(FatalError);
    if (h7.fore()    != (v6 && v6)) FatalErrorInFunction << "test 32f failed" << abort(FatalError);

    h7 = 2.0 && h7;

    if (h7.left()    != (2.0 && (v1 && v1))) FatalErrorInFunction << "test 33a failed" << abort(FatalError);
    if (h7.right()   != (2.0 && (v2 && v2))) FatalErrorInFunction << "test 33b failed" << abort(FatalError);
    if (h7.bottom()  != (2.0 && (v3 && v3))) FatalErrorInFunction << "test 33c failed" << abort(FatalError);
    if (h7.top()     != (2.0 && (v4 && v4))) FatalErrorInFunction << "test 33d failed" << abort(FatalError);
    if (h7.aft()     != (2.0 && (v5 && v5))) FatalErrorInFunction << "test 33e failed" << abort(FatalError);
    if (h7.fore()    != (2.0 && (v6 && v6))) FatalErrorInFunction << "test 33f failed" << abort(FatalError);

    h7 = h7 && 2.0;

    if (h7.left()    != ((2.0 && (v1 && v1)) && 2)) FatalErrorInFunction << "test 34a failed" << abort(FatalError);
    if (h7.right()   != ((2.0 && (v2 && v2)) && 2)) FatalErrorInFunction << "test 34b failed" << abort(FatalError);
    if (h7.bottom()  != ((2.0 && (v3 && v3)) && 2)) FatalErrorInFunction << "test 34c failed" << abort(FatalError);
    if (h7.top()     != ((2.0 && (v4 && v4)) && 2)) FatalErrorInFunction << "test 34d failed" << abort(FatalError);
    if (h7.aft()     != ((2.0 && (v5 && v5)) && 2)) FatalErrorInFunction << "test 34e failed" << abort(FatalError);
    if (h7.fore()    != ((2.0 && (v6 && v6)) && 2)) FatalErrorInFunction << "test 34f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = magSqr(-h7);

    if (h7.left()    != v1*v1) FatalErrorInFunction << "test 35a failed" << abort(FatalError);
    if (h7.right()   != v2*v2) FatalErrorInFunction << "test 35b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3) FatalErrorInFunction << "test 35c failed" << abort(FatalError);
    if (h7.top()     != v4*v4) FatalErrorInFunction << "test 35d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5) FatalErrorInFunction << "test 35e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6) FatalErrorInFunction << "test 35f failed" << abort(FatalError);

    h7 = mag(-h7);

    if (h7.left()    != v1*v1) FatalErrorInFunction << "test 36a failed" << abort(FatalError);
    if (h7.right()   != v2*v2) FatalErrorInFunction << "test 36b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3) FatalErrorInFunction << "test 36c failed" << abort(FatalError);
    if (h7.top()     != v4*v4) FatalErrorInFunction << "test 36d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5) FatalErrorInFunction << "test 36e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6) FatalErrorInFunction << "test 36f failed" << abort(FatalError);

    h7 = cmptMultiply(h7,h7);

    if (h7.left()    != v1*v1*v1*v1) FatalErrorInFunction << "test 37a failed" << abort(FatalError);
    if (h7.right()   != v2*v2*v2*v2) FatalErrorInFunction << "test 37b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3*v3*v3) FatalErrorInFunction << "test 37c failed" << abort(FatalError);
    if (h7.top()     != v4*v4*v4*v4) FatalErrorInFunction << "test 37d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5*v5*v5) FatalErrorInFunction << "test 37e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6*v6*v6) FatalErrorInFunction << "test 37f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = cmptPow(h7,h7);

    if (h7.left()    != Foam::pow(v1,v1)) FatalErrorInFunction << "test 38a failed" << abort(FatalError);
    if (h7.right()   != Foam::pow(v2,v2)) FatalErrorInFunction << "test 38b failed" << abort(FatalError);
    if (h7.bottom()  != Foam::pow(v3,v3)) FatalErrorInFunction << "test 38c failed" << abort(FatalError);
    if (h7.top()     != Foam::pow(v4,v4)) FatalErrorInFunction << "test 38d failed" << abort(FatalError);
    if (h7.aft()     != Foam::pow(v5,v5)) FatalErrorInFunction << "test 38e failed" << abort(FatalError);
    if (h7.fore()    != Foam::pow(v6,v6)) FatalErrorInFunction << "test 38f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = cmptDivide(h7,h7);

    if (h7.left()    != 1) FatalErrorInFunction << "test 39a failed" << abort(FatalError);
    if (h7.right()   != 1) FatalErrorInFunction << "test 39b failed" << abort(FatalError);
    if (h7.bottom()  != 1) FatalErrorInFunction << "test 39c failed" << abort(FatalError);
    if (h7.top()     != 1) FatalErrorInFunction << "test 39d failed" << abort(FatalError);
    if (h7.aft()     != 1) FatalErrorInFunction << "test 39e failed" << abort(FatalError);
    if (h7.fore()    != 1) FatalErrorInFunction << "test 39f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    if (cmptMax(h7) != h7) FatalErrorInFunction << "test 40 failed" << abort(FatalError);
    if (cmptMin(h7) != h7) FatalErrorInFunction << "test 41 failed" << abort(FatalError);

    h7 = cmptSqr(h7);

    if (h7.left()    != v1*v1) FatalErrorInFunction << "test 42a failed" << abort(FatalError);
    if (h7.right()   != v2*v2) FatalErrorInFunction << "test 42b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3) FatalErrorInFunction << "test 42c failed" << abort(FatalError);
    if (h7.top()     != v4*v4) FatalErrorInFunction << "test 42d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5) FatalErrorInFunction << "test 42e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6) FatalErrorInFunction << "test 42f failed" << abort(FatalError);

    h7 = cmptMag(-h7);

    if (h7.left()    != v1*v1) FatalErrorInFunction << "test 43a failed" << abort(FatalError);
    if (h7.right()   != v2*v2) FatalErrorInFunction << "test 43b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3) FatalErrorInFunction << "test 43c failed" << abort(FatalError);
    if (h7.top()     != v4*v4) FatalErrorInFunction << "test 43d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5) FatalErrorInFunction << "test 43e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6) FatalErrorInFunction << "test 43f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = max(h7*h7,h7);

    if (h7.left()    != v1*v1) FatalErrorInFunction << "test 44a failed" << abort(FatalError);
    if (h7.right()   != v2*v2) FatalErrorInFunction << "test 44b failed" << abort(FatalError);
    if (h7.bottom()  != v3*v3) FatalErrorInFunction << "test 44c failed" << abort(FatalError);
    if (h7.top()     != v4*v4) FatalErrorInFunction << "test 44d failed" << abort(FatalError);
    if (h7.aft()     != v5*v5) FatalErrorInFunction << "test 44e failed" << abort(FatalError);
    if (h7.fore()    != v6*v6) FatalErrorInFunction << "test 44f failed" << abort(FatalError);

    h7 = faceScalar(v1,v2,v3,v4,v5,v6);

    h7 = min(h7*h7,h7);

    if (h7.left()    != v1) FatalErrorInFunction << "test 45a failed" << abort(FatalError);
    if (h7.right()   != v2) FatalErrorInFunction << "test 45b failed" << abort(FatalError);
    if (h7.bottom()  != v3) FatalErrorInFunction << "test 45c failed" << abort(FatalError);
    if (h7.top()     != v4) FatalErrorInFunction << "test 45d failed" << abort(FatalError);
    if (h7.aft()     != v5) FatalErrorInFunction << "test 45e failed" << abort(FatalError);
    if (h7.fore()    != v6) FatalErrorInFunction << "test 45f failed" << abort(FatalError);

    minMod(h3,h4);

    if (h7.lower().x() != v1) FatalErrorInFunction << "test 46a failed" << abort(FatalError);
    if (h7.lower().y() != v3) FatalErrorInFunction << "test 46b failed" << abort(FatalError);
    if (h7.lower().z() != v5) FatalErrorInFunction << "test 46c failed" << abort(FatalError);

    if (h7.upper().x() != v2) FatalErrorInFunction << "test 46a failed" << abort(FatalError);
    if (h7.upper().y() != v4) FatalErrorInFunction << "test 46b failed" << abort(FatalError);
    if (h7.upper().z() != v6) FatalErrorInFunction << "test 46c failed" << abort(FatalError);

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
    asin(h3);
    acos(h3);
    atan(h3);
    sinh(h3);
    cosh(h3);
    tanh(h3);
    asinh(h3);
    acosh(h3);
    atanh(h3);
    erf(h3);
    erfc(h3);
    lgamma(h3);
    j0(h3);
    j1(h3);
    y0(h3);
    y1(h3);
}
