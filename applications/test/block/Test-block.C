#include "arguments.H"

#include "block.H"

#include "fileOperation.H"
#include "OSspecific.H"

#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;
using namespace briscola;

template<class Type>
void testConstructors()
{
    block<Type> b1(2,3,4);
    block<Type> b2(2,3,4,Zero);
    block<Type> b3(2,3,4,pTraits<Type>::one);

    block<Type> b4(labelVector(2,3,4));
    block<Type> b5(labelVector(2,3,4),Zero);
    block<Type> b6(labelVector(2,3,4),pTraits<Type>::one);

    List<Type> list(24, pTraits<Type>::one*2);

    block<Type> b7(2,3,4,list);
    block<Type> b8(labelVector(2,3,4),list);

    Type arr[24];

    for (label i = 0; i < 24; i++)
    {
        arr[i] = pTraits<Type>::one*i;
    }

    block<Type> b9(2,3,4, arr);
    block<Type> b10(labelVector(2,3,4), arr);

    block<Type> b11a(b10);
    block<Type> b11b(b10, Zero);
    block<Type> b11c(b10, pTraits<Type>::one*2);

    block<Type> b12a(b10*2.0);
    block<Type> b12b(b10*2.0, Zero);
    block<Type> b12c(b10*2.0, pTraits<Type>::one*2);

    const word fileName =
        "dummy-"
      + word(pTraits<Type>::typeName)
      + "-"
      + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << b10 << endl;

    IFstream is(fileName);

    rm(fileName);

    block<Type> b13(is);

    block<Type> b14;

    b14.transferData(b13);
}

template<class Type>
void testIndexing()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
    }

    if (b1.l() != shape.x()) FatalErrorInFunction << "test 1 failed" << abort(FatalError);
    if (b1.m() != shape.y()) FatalErrorInFunction << "test 2 failed" << abort(FatalError);
    if (b1.n() != shape.z()) FatalErrorInFunction << "test 3 failed" << abort(FatalError);

    if (b1.cbegin() != &b1(0))
        FatalErrorInFunction << "test 4 failed" << abort(FatalError);

    if (b1.size() != cmptProduct(shape))
        FatalErrorInFunction << "test 5 failed" << abort(FatalError);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 6 failed" << abort(FatalError);

    if (b1.shapeTensor() != labelTensor(shape.x(),0,0,0,shape.y(),0,0,0,shape.z()))
        FatalErrorInFunction << "test 7 failed" << abort(FatalError);

    l = 0;

    forAllBlockLinear(b1, i)
    {
        if (b1(i) != pTraits<Type>::one*l++)
            FatalErrorInFunction << "test 8 failed" << abort(FatalError);
    }

    l = 0;

    forAllBlock(b1, i, j, k)
    {
        if (b1(i,j,k) != pTraits<Type>::one*l++)
            FatalErrorInFunction << "test 8 failed" << abort(FatalError);
    }

    // Shifted interpolation

    for (int i = 0; i < shape.x()-1; i++)
    for (int j = 1; j < shape.y(); j++)
    for (int k = 0; k < shape.z()-1; k++)
    {
        b1.interp(i,j,k,vector(0.1,-0.2,0.3));
        b1.interp(labelVector(i,j,k),vector(0.1,-0.2,0.3));
    }
}

template<class Type>
void testTransformations()
{
    labelVector shape(2,3,4);

    block<Type> b1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
    }

    // Reflections shouldn't change shape

    b1.transform(reflectX);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    b1.transform(reflectY);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    b1.transform(reflectZ);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    // 180 degree rotations shouldn't change shape

    b1.transform(rotateX2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    b1.transform(rotateY2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    b1.transform(rotateZ2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    // Permutations

    b1.transform(permuteXY);
    shape = permuteXY & shape;

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    b1.transform(permuteXZ);
    shape = permuteXZ & shape;

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    b1.transform(permuteYZ);
    shape = permuteYZ & shape;

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    // Rotations

    b1.transform(rotateX1);
    shape = cmptMag(rotateX1 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4a failed" << abort(FatalError);

    b1.transform(rotateX3);
    shape = cmptMag(rotateX3 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4b failed" << abort(FatalError);

    b1.transform(rotateY1);
    shape = cmptMag(rotateY1 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4c failed" << abort(FatalError);

    b1.transform(rotateY3);
    shape = cmptMag(rotateY3 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4d failed" << abort(FatalError);

    b1.transform(rotateZ1);
    shape = cmptMag(rotateZ1 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4e failed" << abort(FatalError);

    b1.transform(rotateZ3);
    shape = cmptMag(rotateZ3 & shape);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 4f failed" << abort(FatalError);

    // Slicing

    block<Type> b2(b1.slice(1,1));
    shape.y() = 1;

    if (b2.shape() != shape)
        FatalErrorInFunction << "test 5 failed" << abort(FatalError);

    b2.squeeze();
    shape = labelVector(shape.x(), shape.z(), 1);

    if (b2.shape() != shape)
        FatalErrorInFunction << "test 6 failed" << abort(FatalError);
}

template<class Type>
void testMemberOperators()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);
    block<scalar> s1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        s1(i,j,k) = scalar(2.0*l+1.0);
        b1(i,j,k) = Zero;
        b2(i,j,k) = pTraits<Type>::one*l++;
    }

    b1 = Zero;
    b1 = pTraits<Type>::one*2.0;
    b1 = b2;
    b1 = (2*b2);

    b1 += b2;
    b1 += (2*b2);

    b1 -= b2;
    b1 -= (2*b2);

    b1 *= s1;
    b1 *= (2*s1);

    b1 /= s1;
    b1 /= (2*s1);

    b1 += pTraits<Type>::one;
    b1 -= pTraits<Type>::one;

    b1 *= scalar(2.0);
    b1 /= scalar(2.0);
}

template<class Type>
void testPrimitiveFunctions()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);
    block<scalar> s1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
        b2(i,j,k) = pTraits<Type>::one*l++;
        s1(i,j,k) = (l+++1);
    }

    mag(b1);
    mag(b1*2.0);

    max(b1);
    max(b1*2.0);
    gMax(b1);
    gMax(b1*2.0);

    min(b1);
    min(b1*2.0);
    gMin(b1);
    gMin(b1*2.0);

    sum(b1);
    sum(b1*2.0);
    gSum(b1);
    gSum(b1*2.0);

    average(b1);
    average(b1*2.0);
    gAverage(b1);
    gAverage(b1*2.0);

    sumProd(b1, b1);
    sumProd(b1*2.0, b1);
    sumProd(b1, b1*2.0);
    sumProd(b1*2.0, b1*2.0);

    max(b1,b2);
    max(b1*2.0,b2);
    max(b1,b2*2.0);
    max(b1*2.0,b2*2.0);

    min(b1,b2);
    min(b1*2.0,b2);
    min(b1,b2*2.0);
    min(b1*2.0,b2*2.0);

    max(b1,pTraits<Type>::one);
    max(b1*2.0,pTraits<Type>::one);
    max(pTraits<Type>::one,b2);
    max(pTraits<Type>::one,b2*2.0);

    min(b1,pTraits<Type>::one);
    min(b1*2.0,pTraits<Type>::one);
    min(pTraits<Type>::one,b2);
    min(pTraits<Type>::one,b2*2.0);

    -b2;

    b1*s1;
    (b1*2.0)*s1;
    s1*b1;
    s1*(b1*2.0);

    b1/s1;
    (b1*2.0)/s1;

    b1+b2;

    b1-b2;
}

template<class Type>
void testVectorSpaceFunctions()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
        b2(i,j,k) = pTraits<Type>::one*(l+++1);
    }

    cmptMax(b1);
    cmptMax(b1*2.0);

    cmptMin(b1);
    cmptMin(b1*2.0);

    cmptAv(b1);
    cmptAv(b1*2.0);

    cmptMag(b1);
    cmptMag(b1*2.0);

    maxMagSqr(b1);
    maxMagSqr(b1*2.0);
    gMaxMagSqr(b1);
    gMaxMagSqr(b1*2.0);

    minMagSqr(b1);
    minMagSqr(b1*2.0);
    gMinMagSqr(b1);
    gMinMagSqr(b1*2.0);

    sumMag(b1);
    sumMag(b1*2.0);
    gSumMag(b1);
    gSumMag(b1*2.0);

    sumCmptProd(b1, b1);
    sumCmptProd(b1*2.0, b1);
    sumCmptProd(b1, b1*2.0);
    sumCmptProd(b1*2.0, b1*2.0);

    sumCmptMag(b1);
    sumCmptMag(b1*2.0);
    gSumCmptMag(b1);
    gSumCmptMag(b1*2.0);

    cmptMultiply(b1,b2);
    cmptMultiply(b1*2.0,b2);
    cmptMultiply(b1,b2*2.0);
    cmptMultiply(b1*2.0,b2*2.0);

    cmptDivide(b1,b2);
    cmptDivide(b1*2.0,b2);
    cmptDivide(b1,b2*2.0);
    cmptDivide(b1*2.0,b2*2.0);

    cmptMultiply(b1,pTraits<Type>::one);
    cmptMultiply(b1*2.0,pTraits<Type>::one);
    cmptMultiply(pTraits<Type>::one,b2);
    cmptMultiply(pTraits<Type>::one,b2*2.0);

    cmptDivide(b1,pTraits<Type>::one);
    cmptDivide(b1*2.0,pTraits<Type>::one);
    cmptDivide(pTraits<Type>::one,b2);
    cmptDivide(pTraits<Type>::one,b2*2.0);
}

template<class Type>
void testStencilFunctions()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);
    block<scalar> s1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
        b2(i,j,k) = pTraits<Type>::one*(l+++1);
        s1(i,j,k) = (l+++1);
    }

    -b1;

    b1+b2;
    (2.0*b1)+b2;
    b1+(2.0*b2);
    (2.0*b1)+(2.0*b2);

    b1-b2;
    (2.0*b1)-b2;
    b1-(2.0*b2);
    (2.0*b1)-(2.0*b2);

    b1-b2;
    (2.0*b1)-b2;
    b1-(2.0*b2);
    (2.0*b1)-(2.0*b2);

    b1*s1;
    s1*b1;
    (2.0*b1)*s1;
    s1*(2.0*b1);

    b1*(2.0*s1);
    (2.0*s1)*b1;
    (2.0*b1)*(2.0*s1);
    (2.0*s1)*(2.0*b1);

    b1/s1;
    (2.0*b1)/s1;

    b1/(2.0*s1);
    (2.0*b1)/(2.0*s1);
}

void testScalarFunctions()
{
    const labelVector shape(2,3,4);

    block<scalar> b1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = scalar(l+++1);
    }

    sumSqr(b1);
    gSumSqr(b1);

    b1/b1;
    (b1*2.0)/b1;
    b1/(b1*2.0);
}

void testVectorFunctions()
{
    const labelVector shape(2,3,4);

    block<vector> b1(shape);
    block<vector> b2(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = vector(l,l+1,l+2); l++;
        b2(i,j,k) = vector(l,l+1,l+2); l++;
    }

    sumSqr(b1);
    sumSqr(b1*2.0);
    gSumSqr(b1);
    gSumSqr(b1*2.0);

    b1*b2;
    (b1*2.0)*b2;
    b1*(b2*2.0);
    (b1*2.0)*(b2*2.0);

    b1 & b2;
    (b1*2.0) & b2;
    b1 & (b2*2.0);
    (b1*2.0) & (b2*2.0);

    b1 ^ b2;
    (b1*2.0) ^ b2;
    b1 ^ (b2*2.0);
    (b1*2.0) ^ (b2*2.0);
}

void testTensorFunctions()
{
    const labelVector shape(2,3,4);

    block<tensor> b1(shape);
    block<symmTensor> b2(shape);
    block<sphericalTensor> b3(shape);
    block<diagTensor> b4(shape);
    block<vector> v1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<tensor>::one*l++;
        b2(i,j,k) = pTraits<symmTensor>::one*l++;
        b3(i,j,k) = pTraits<sphericalTensor>::one*l++;
        b4(i,j,k) = pTraits<diagTensor>::one*l++;
        v1(i,j,k) = vector(l,l+1,l+2); l++;
    }

    b1 & b1;
    (b1*2.0) & b1;
    b1 & (b1*2.0);
    (b1*2.0) & (b1*2.0);

    b1 & b2;
    (b1*2.0) & b2;
    b1 & (b2*2.0);
    (b1*2.0) & (b2*2.0);

    b1 && b1;
    (b1*2.0) && b1;
    b1 && (b1*2.0);
    (b1*2.0) && (b1*2.0);

    b1 && b2;
    (b1*2.0) && b2;
    b1 && (b2*2.0);
    (b1*2.0) && (b2*2.0);

    b1 & v1;
    (b1*2.0) & v1;
    v1 & (b1*2.0);
    (v1*2.0) & (b1*2.0);

    b2 & v1;
    (b2*2.0) & v1;
    v1 & (b2*2.0);
    (v1*2.0) & (b2*2.0);

    b3 & v1;
    (b3*2.0) & v1;
    v1 & (b3*2.0);
    (v1*2.0) & (b3*2.0);

    b4 & v1;
    (b4*2.0) & v1;
    v1 & (b4*2.0);
    (v1*2.0) & (b4*2.0);
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    testConstructors<label>();
    testConstructors<scalar>();
    testConstructors<hexScalar>();
    testConstructors<vector>();
    testConstructors<hexVector>();
    testConstructors<tensor>();
    testConstructors<symmTensor>();
    testConstructors<sphericalTensor>();
    testConstructors<diagTensor>();
    testConstructors<stencil>();
    testConstructors<symmStencil>();
    testConstructors<diagStencil>();

    testIndexing<label>();
    testIndexing<scalar>();
    testIndexing<hexScalar>();
    testIndexing<vector>();
    testIndexing<hexVector>();
    testIndexing<tensor>();
    testIndexing<symmTensor>();
    testIndexing<sphericalTensor>();
    testIndexing<diagTensor>();
    testIndexing<stencil>();
    testIndexing<symmStencil>();
    testIndexing<diagStencil>();

    testTransformations<label>();
    testTransformations<scalar>();
    testTransformations<hexScalar>();
    testTransformations<vector>();
    testTransformations<hexVector>();
    testTransformations<tensor>();
    testTransformations<symmTensor>();
    testTransformations<sphericalTensor>();
    testTransformations<diagTensor>();
    testTransformations<stencil>();
    testTransformations<symmStencil>();
    testTransformations<diagStencil>();

    testMemberOperators<label>();
    testMemberOperators<scalar>();
    testMemberOperators<hexScalar>();
    testMemberOperators<vector>();
    testMemberOperators<hexVector>();
    testMemberOperators<tensor>();
    testMemberOperators<symmTensor>();
    testMemberOperators<sphericalTensor>();
    testMemberOperators<diagTensor>();
    testMemberOperators<stencil>();
    testMemberOperators<symmStencil>();
    testMemberOperators<diagStencil>();

    testPrimitiveFunctions<scalar>();
    testPrimitiveFunctions<vector>();
    testPrimitiveFunctions<tensor>();
    testPrimitiveFunctions<symmTensor>();
    testPrimitiveFunctions<sphericalTensor>();
    testPrimitiveFunctions<diagTensor>();

    testVectorSpaceFunctions<vector>();
    testVectorSpaceFunctions<tensor>();
    testVectorSpaceFunctions<symmTensor>();
    testVectorSpaceFunctions<sphericalTensor>();
    testVectorSpaceFunctions<diagTensor>();

    testStencilFunctions<stencil>();
    testStencilFunctions<symmStencil>();
    testStencilFunctions<diagStencil>();

    // Type-specific functions

    testScalarFunctions();
    testVectorFunctions();
    testTensorFunctions();
}
