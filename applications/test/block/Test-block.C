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

    block<Type> b13a(false, b10);
    block<Type> b13b(false, b10, Zero);
    block<Type> b13c(false, b10, pTraits<Type>::one*2);

    block<Type> b14a(true, b13a);
    block<Type> b14b(true, b13b, Zero);
    block<Type> b14c(true, b13c, pTraits<Type>::one*2);

    const word fileName =
        "dummy-"
      + word(pTraits<Type>::typeName)
      + "-"
      + Foam::name(Pstream::myProcNo());

    OFstream os(fileName);

    os << b10 << endl;

    IFstream is(fileName);

    rm(fileName);

    block<Type> b14(is);

    block<Type> b15;

    b14.transfer(b15);
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

    // Interpolation to a scalar index

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = (i+2*j+3*k)*pTraits<Type>::one;
    }

    for (int i = 0; i < shape.x()-1; i++)
    for (int j = 1; j < shape.y(); j++)
    for (int k = 0; k < shape.z()-1; k++)
    {
        const scalar ii = i+0.1;
        const scalar jj = j-0.2;
        const scalar kk = k+0.3;

        Type r1 = b1.interp(ii,jj,kk);
        Type r2 = b1.interp(vector(ii,jj,kk));

        if ((r1-(ii+2*jj+3*kk)*pTraits<Type>::one) > 1e-12*pTraits<Type>::one)
        {
            FatalErrorInFunction << "test 9a failed" << abort(FatalError);
        }

        if ((r2-(ii+2*jj+3*kk)*pTraits<Type>::one) > 1e-12*pTraits<Type>::one)
        {
            FatalErrorInFunction << "test 9b failed" << abort(FatalError);
        }
    }
}

template<class Type>
void testTransformations()
{
    const label l = 2;
    const label m = 3;
    const label n = 4;

    const labelVector shape(l,m,n);

    block<Type> b1(shape);

    label c = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*c++;
    }

    const block<Type> b0(b1);

    // Reflections shouldn't change shape

    b1 = b0;
    b1.transform(reflectX);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(l-1-i,j,k))
            FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    b1 = b0;
    b1.transform(reflectY);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,m-1-j,k))
            FatalErrorInFunction << "test 1d failed" << abort(FatalError);

    b1 = b0;
    b1.transform(reflectZ);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 1e failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,n-1-k))
            FatalErrorInFunction << "test 1f failed" << abort(FatalError);

    // 180 degree rotations shouldn't change shape

    b1 = b0;
    b1.transform(rotateX2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,m-1-j,n-1-k))
            FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateY2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(l-1-i,j,n-1-k))
            FatalErrorInFunction << "test 2d failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateZ2);

    if (b1.shape() != shape)
        FatalErrorInFunction << "test 2e failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(l-1-i,m-1-j,k))
            FatalErrorInFunction << "test 2f failed" << abort(FatalError);

    // Permutations

    b1 = b0;
    b1.transform(permuteXY);

    if (b1.shape() != (permuteXY & shape))
        FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(j,i,k))
            FatalErrorInFunction << "test 3b failed" << abort(FatalError);

    b1 = b0;
    b1.transform(permuteXZ);

    if (b1.shape() != (permuteXZ & shape))
        FatalErrorInFunction << "test 3c failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(k,j,i))
            FatalErrorInFunction << "test 3d failed" << abort(FatalError);

    b1 = b0;
    b1.transform(permuteYZ);

    if (b1.shape() != (permuteYZ & shape))
        FatalErrorInFunction << "test 3e failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,k,j))
            FatalErrorInFunction << "test 3f failed" << abort(FatalError);

    // Rotations

    b1 = b0;
    b1.transform(rotateX1);

    if (b1.shape() != cmptMag(rotateX1 & shape))
        FatalErrorInFunction << "test 4a failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,k,n-1-j))
            FatalErrorInFunction << "test 4b failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateX3);

    if (b1.shape() != cmptMag(rotateX3 & shape))
        FatalErrorInFunction << "test 4c failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,m-1-k,j))
            FatalErrorInFunction << "test 4d failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateY1);

    if (b1.shape() != cmptMag(rotateY1 & shape))
        FatalErrorInFunction << "test 4e failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(l-1-k,j,i))
            FatalErrorInFunction << "test 4f failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateY3);

    if (b1.shape() != cmptMag(rotateY3 & shape))
        FatalErrorInFunction << "test 4g failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(k,j,n-1-i))
            FatalErrorInFunction << "test 4h failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateZ1);

    if (b1.shape() != cmptMag(rotateZ1 & shape))
        FatalErrorInFunction << "test 4i failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(j,m-1-i,k))
            FatalErrorInFunction << "test 4j failed" << abort(FatalError);

    b1 = b0;
    b1.transform(rotateZ3);

    if (b1.shape() != cmptMag(rotateZ3 & shape))
        FatalErrorInFunction << "test 4k failed" << abort(FatalError);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(l-1-j,i,k))
            FatalErrorInFunction << "test 4l failed" << abort(FatalError);

    // Slicing

    b1 = b0;

    for (label dir = 0; dir < 3; dir++)
    {
        block<Type> s(b1.slice(dir+1,dir));
        labelVector shape2(shape);
        shape2[dir] = 1;

        if (s.shape() != shape2)
            FatalErrorInFunction << "test 5 failed" << abort(FatalError);

        forAllBlock(s, i, j, k)
        {
            const labelVector ijk = units[dir]*(1+dir) + labelVector(i,j,k);

            if (s(i,j,k) != b0(ijk))
                FatalErrorInFunction << "test 6 failed" << abort(FatalError);
        }
    }
}

template<class Type>
void testMemberOperators()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<scalar> s1(shape);

    label l = 0;

    forAllBlock(b1, i, j, k)
    {
        s1(i,j,k) = scalar(2.0*l+1.0);
        b1(i,j,k) = pTraits<Type>::one*l++;
    }

    const block<Type> b0(b1);

    b1 = Zero;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != pTraits<Type>::zero)
            FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    b1 = pTraits<Type>::one*2.0;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != pTraits<Type>::one*2.0)
            FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    b1 = b0;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k))
            FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    b1 = (2*b0);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != 2*b0(i,j,k))
            FatalErrorInFunction << "test 1d failed" << abort(FatalError);

    b1 = b0;
    b1 += b0;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != 2*b0(i,j,k))
            FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    b1 += (2*b0);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != 4*b0(i,j,k))
            FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    b1 -= b0;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != 3*b0(i,j,k))
            FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    b1 -= (2*b0);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k))
            FatalErrorInFunction << "test 2d failed" << abort(FatalError);

    b1 = b0;
    b1 *= s1;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    b1 *= (2*s1);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k)*s1(i,j,k)*s1(i,j,k)*2)
            FatalErrorInFunction << "test 3b failed" << abort(FatalError);

    1.0/s1;
    1.0/(1.0*s1);

    b1 /= s1;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k)*s1(i,j,k)*2)
            FatalErrorInFunction << "test 3c failed" << abort(FatalError);

    b1 /= (2*s1);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k))
            FatalErrorInFunction << "test 3d failed" << abort(FatalError);

    b1 = b0;

    b1 += pTraits<Type>::one;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k) + pTraits<Type>::one)
            FatalErrorInFunction << "test 4a failed" << abort(FatalError);

    b1 -= pTraits<Type>::one;

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k))
            FatalErrorInFunction << "test 4b failed" << abort(FatalError);

    b1 *= scalar(2.0);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != 2*b0(i,j,k))
            FatalErrorInFunction << "test 5a failed" << abort(FatalError);

    b1 /= scalar(2.0);

    forAllBlock(b1, i, j, k)
        if (b1(i,j,k) != b0(i,j,k))
            FatalErrorInFunction << "test 5b failed" << abort(FatalError);
}

template<class Type>
void testPrimitiveFunctions()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);
    block<scalar> s1(shape);

    label l = 3;

    Type sumb1(Zero);

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = pTraits<Type>::one*l++;
        b2(i,j,k) = pTraits<Type>::one*l++;
        s1(i,j,k) = (l+++1);

        sumb1 += b1(i,j,k);
    }

    const Type avb1 = sumb1/cmptProduct(shape);

    mag(b1);
    mag(b1*2.0);

    if (max(b1) != pTraits<Type>::one*cmptProduct(shape)*3)
        FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    if (max(b1*2.0) != pTraits<Type>::one*cmptProduct(shape)*3*2)
        FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    if (gMax(b1) != pTraits<Type>::one*cmptProduct(shape)*3)
        FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    if (gMax(b1*2.0) != pTraits<Type>::one*cmptProduct(shape)*3*2)
        FatalErrorInFunction << "test 1d failed" << abort(FatalError);


    if (min(b1) != pTraits<Type>::one*3)
        FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    if (min(b1*2.0) != pTraits<Type>::one*3*2)
        FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    if (gMin(b1) != pTraits<Type>::one*3)
        FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    if (gMin(b1*2.0) != pTraits<Type>::one*3*2)
        FatalErrorInFunction << "test 2d failed" << abort(FatalError);


    if (sum(b1) != sumb1)
        FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    if (sum(b1*2.0) != sumb1*2)
        FatalErrorInFunction << "test 3b failed" << abort(FatalError);

    if (gSum(b1) != sumb1)
        FatalErrorInFunction << "test 3c failed" << abort(FatalError);

    if (gSum(b1*2) != sumb1*2)
        FatalErrorInFunction << "test 3d failed" << abort(FatalError);


    if (average(b1) != avb1)
        FatalErrorInFunction << "test 4a failed" << abort(FatalError);

    if (average(b1*2) != avb1*2)
        FatalErrorInFunction << "test 4b failed" << abort(FatalError);

    if (gAverage(b1) != avb1)
        FatalErrorInFunction << "test 4c failed" << abort(FatalError);

    if (gAverage(b1*2) != avb1*2)
        FatalErrorInFunction << "test 4d failed" << abort(FatalError);


    sumProd(b1, b1);
    sumProd(b1*2.0, b1);
    sumProd(b1, b1*2.0);
    sumProd(b1*2.0, b1*2.0);


    block<Type> b3 = max(b1,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b2(i,j,k))
            FatalErrorInFunction << "test 5a failed" << abort(FatalError);

    b3 = max(b1*2.0,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k))
            FatalErrorInFunction << "test 5b failed" << abort(FatalError);

    b3 = max(b1,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b2(i,j,k))
            FatalErrorInFunction << "test 5c failed" << abort(FatalError);

    b3 = max(b1*2.0,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b2(i,j,k))
            FatalErrorInFunction << "test 5d failed" << abort(FatalError);


    b3 = min(b1,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 6a failed" << abort(FatalError);

    b3 = min(b1*2.0,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b2(i,j,k))
            FatalErrorInFunction << "test 6b failed" << abort(FatalError);

    b3 = min(b1,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 6c failed" << abort(FatalError);

    b3 = min(b1*2.0,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k))
            FatalErrorInFunction << "test 6d failed" << abort(FatalError);


    b3 = max(b1,pTraits<Type>::one);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 7a failed" << abort(FatalError);

    b3 = max(b1*2.0,pTraits<Type>::one);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k))
            FatalErrorInFunction << "test 7b failed" << abort(FatalError);

    b3 = max(pTraits<Type>::one,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != b2(i,j,k))
            FatalErrorInFunction << "test 7c failed" << abort(FatalError);

    b3 = max(pTraits<Type>::one,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != 2*b2(i,j,k))
            FatalErrorInFunction << "test 7d failed" << abort(FatalError);


    b3 = min(b1,pTraits<Type>::one);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one)
            FatalErrorInFunction << "test 8a failed" << abort(FatalError);

    b3 = min(b1*2.0,pTraits<Type>::one);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one)
            FatalErrorInFunction << "test 8b failed" << abort(FatalError);

    b3 = min(pTraits<Type>::one,b2);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one)
            FatalErrorInFunction << "test 8c failed" << abort(FatalError);

    b3 = min(pTraits<Type>::one,b2*2.0);

    forAllBlock(b1, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one)
            FatalErrorInFunction << "test 8d failed" << abort(FatalError);


    -b2;

    b2 = b1*s1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9a failed" << abort(FatalError);

    b2 = (b1*2.0)*s1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9b failed" << abort(FatalError);

    b2 = b1*(2.0*s1);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9c failed" << abort(FatalError);

    b2 = s1*b1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9d failed" << abort(FatalError);

    b2 = s1*(b1*2.0);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9e failed" << abort(FatalError);

    b2 = (2.0*s1)*b1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k)*s1(i,j,k))
            FatalErrorInFunction << "test 9f failed" << abort(FatalError);


    b2 = b1/s1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k)/s1(i,j,k))
            FatalErrorInFunction << "test 10a failed" << abort(FatalError);

    b2 = (b1*2.0)/s1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k)/s1(i,j,k))
            FatalErrorInFunction << "test 10b failed" << abort(FatalError);

    b2 = b1/(2.0*s1);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 0.5*b1(i,j,k)/s1(i,j,k))
            FatalErrorInFunction << "test 10c failed" << abort(FatalError);


    b3 = b1+b2;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)+b2(i,j,k))
            FatalErrorInFunction << "test 11a failed" << abort(FatalError);

    b3 = (2.0*b1)+b2;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k)+b2(i,j,k))
            FatalErrorInFunction << "test 11b failed" << abort(FatalError);

    b3 = b1+(2.0*b2);

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)+2*b2(i,j,k))
            FatalErrorInFunction << "test 11c failed" << abort(FatalError);

    b3 = b1-b2;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)-b2(i,j,k))
            FatalErrorInFunction << "test 12a failed" << abort(FatalError);

    b3 = (2.0*b1)-b2;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k)-b2(i,j,k))
            FatalErrorInFunction << "test 12b failed" << abort(FatalError);

    b3 = b1-(2.0*b2);

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)-2*b2(i,j,k))
            FatalErrorInFunction << "test 12c failed" << abort(FatalError);


    b3 = b1 + pTraits<Type>::one;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)+pTraits<Type>::one)
            FatalErrorInFunction << "test 13a failed" << abort(FatalError);

    b3 = 2*b1 + pTraits<Type>::one;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k)+pTraits<Type>::one)
            FatalErrorInFunction << "test 13b failed" << abort(FatalError);

    b3 = b1 - pTraits<Type>::one;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)-pTraits<Type>::one)
            FatalErrorInFunction << "test 13c failed" << abort(FatalError);

    b3 = 2*b1 - pTraits<Type>::one;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k)-pTraits<Type>::one)
            FatalErrorInFunction << "test 13d failed" << abort(FatalError);

    b3 = pTraits<Type>::one + b1;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != b1(i,j,k)+pTraits<Type>::one)
            FatalErrorInFunction << "test 13e failed" << abort(FatalError);

    b3 = pTraits<Type>::one + 2*b1;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != 2*b1(i,j,k)+pTraits<Type>::one)
            FatalErrorInFunction << "test 13f failed" << abort(FatalError);

    b3 = pTraits<Type>::one - b1;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one - b1(i,j,k))
            FatalErrorInFunction << "test 13g failed" << abort(FatalError);

    b3 = pTraits<Type>::one - 2*b1;

    forAllBlock(b2, i, j, k)
        if (b3(i,j,k) != pTraits<Type>::one - 2*b1(i,j,k))
            FatalErrorInFunction << "test 13h failed" << abort(FatalError);
}

template<class Type>
void testVectorSpaceFunctions()
{
    const labelVector shape(2,3,4);

    block<Type> b1(shape);
    block<Type> b2(shape);

    const label nComp = pTraits<Type>::nComponents;

    scalar b1sm = 0.0;
    Type b1b1scp(Zero);
    Type b1scm(Zero);

    forAllBlock(b1, i, j, k)
    {
        Type a;

        for (label dir = 0; dir < nComp; dir++)
        {
            a[dir] = i*j*k+dir;
        }

        b1(i,j,k) = a;
        b2(i,j,k) = a*2;

        b1sm += Foam::mag(a);
        b1b1scp += Foam::cmptMultiply(a,a);
        b1scm += Foam::cmptMag(a);
    }


    scalarBlock r1 = cmptMax(b1);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != i*j*k+(nComp-1))
            FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    r1 = cmptMax(b1*2.0);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != 2*(i*j*k+(nComp-1)))
            FatalErrorInFunction << "test 1b failed" << abort(FatalError);


    r1 = cmptMin(b1);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != i*j*k)
            FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    r1 = cmptMin(b1*2.0);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != 2*i*j*k)
            FatalErrorInFunction << "test 2b failed" << abort(FatalError);


    r1 = cmptAv(b1);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != i*j*k + (nComp-1.0)/2.0)
            FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    r1 = cmptAv(b1*2.0);

    forAllBlock(r1, i, j, k)
        if (r1(i,j,k) != 2.0*(i*j*k + (nComp-1.0)/2.0))
            FatalErrorInFunction << "test 3b failed" << abort(FatalError);


    r1 = mag(b1);

    forAllBlock(r1, i, j, k)
    {
        if (r1(i,j,k) != Foam::mag(b1(i,j,k)))
            FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    }

    r1 = mag(2*b1);

    forAllBlock(r1, i, j, k)
    {
        if (r1(i,j,k) != Foam::mag(2.0*b1(i,j,k)))
            FatalErrorInFunction << "test 4a failed" << abort(FatalError);
    }


    b2 = cmptMag(-b1);

    forAllBlock(r1, i, j, k)
        if (b2(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 5a failed" << abort(FatalError);

    b2 = cmptMag(-b1*2.0);

    forAllBlock(r1, i, j, k)
        if (b2(i,j,k) != 2*b1(i,j,k))
            FatalErrorInFunction << "test 5b failed" << abort(FatalError);


    if (maxMagSqr(b1) != Foam::magSqr(b1(b1.size()-1)))
        FatalErrorInFunction << "test 6a failed" << abort(FatalError);

    if (maxMagSqr(2.0*b1) != Foam::magSqr(2.0*b1(b1.size()-1)))
        FatalErrorInFunction << "test 6b failed" << abort(FatalError);

    if (gMaxMagSqr(b1) != Foam::magSqr(b1(b1.size()-1)))
        FatalErrorInFunction << "test 6c failed" << abort(FatalError);

    if (gMaxMagSqr(2.0*b1) != Foam::magSqr(2.0*b1(b1.size()-1)))
        FatalErrorInFunction << "test 6d failed" << abort(FatalError);


    if (minMagSqr(b1) != Foam::magSqr(b1(0)))
        FatalErrorInFunction << "test 7a failed" << abort(FatalError);

    if (minMagSqr(2.0*b1) != Foam::magSqr(2.0*b1(0)))
        FatalErrorInFunction << "test 7b failed" << abort(FatalError);

    if (gMinMagSqr(b1) != Foam::magSqr(b1(0)))
        FatalErrorInFunction << "test 7c failed" << abort(FatalError);

    if (gMinMagSqr(2.0*b1) != Foam::magSqr(2.0*b1(0)))
        FatalErrorInFunction << "test 7d failed" << abort(FatalError);


    if (sumMag(b1) != b1sm)
        FatalErrorInFunction << "test 8a failed" << abort(FatalError);

    if (sumMag(1.0*b1) != b1sm)
        FatalErrorInFunction << "test 8b failed" << abort(FatalError);

    if (gSumMag(b1) != b1sm)
        FatalErrorInFunction << "test 8c failed" << abort(FatalError);

    if (gSumMag(1.0*b1) != b1sm)
        FatalErrorInFunction << "test 8d failed" << abort(FatalError);


    if (sumCmptProd(b1, b1) != b1b1scp)
        FatalErrorInFunction << "test 9a failed" << abort(FatalError);

    if (sumCmptProd(2.0*b1, b1) != 2.0*b1b1scp)
        FatalErrorInFunction << "test 9b failed" << abort(FatalError);

    if (sumCmptProd(b1, b1*2.0) != 2.0*b1b1scp)
        FatalErrorInFunction << "test 9c failed" << abort(FatalError);

    if (sumCmptProd(2.0*b1, 2.0*b1) != 4.0*b1b1scp)
        FatalErrorInFunction << "test 9d failed" << abort(FatalError);


    if (sumCmptMag(b1) != b1scm)
        FatalErrorInFunction << "test 10a failed" << abort(FatalError);

    if (sumCmptMag(2.0*b1) != 2.0*b1scm)
        FatalErrorInFunction << "test 10b failed" << abort(FatalError);

    if (gSumCmptMag(b1) != b1scm)
        FatalErrorInFunction << "test 10c failed" << abort(FatalError);

    if (gSumCmptMag(2.0*b1) != 2.0*b1scm)
        FatalErrorInFunction << "test 10d failed" << abort(FatalError);


    b2 = cmptMultiply(b1,b1);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptMultiply(b1(i,j,k),b1(i,j,k)))
            FatalErrorInFunction << "test 11a failed" << abort(FatalError);

    b2 = cmptMultiply(b1*2.0,b1);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2.0*cmptMultiply(b1(i,j,k),b1(i,j,k)))
            FatalErrorInFunction << "test 11b failed" << abort(FatalError);

    b2 = cmptMultiply(b1,b1*2.0);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2.0*cmptMultiply(b1(i,j,k),b1(i,j,k)))
            FatalErrorInFunction << "test 11c failed" << abort(FatalError);

    b2 = cmptMultiply(b1*2.0,b1*2.0);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 4.0*cmptMultiply(b1(i,j,k),b1(i,j,k)))
            FatalErrorInFunction << "test 11d failed" << abort(FatalError);


    b2 = cmptDivide(b1,(b1+pTraits<Type>::one)());

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(b1(i,j,k), b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 12a failed" << abort(FatalError);

    b2 = cmptDivide(b1*2.0,(b1+pTraits<Type>::one)());

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(2.0*b1(i,j,k), b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 12b failed" << abort(FatalError);

    b2 = cmptDivide(b1,b1+pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(b1(i,j,k), b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 12c failed" << abort(FatalError);

    b2 = cmptDivide(b1*2.0,b1+pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(2.0*b1(i,j,k), b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 12d failed" << abort(FatalError);


    b2 = cmptMultiply(b1,pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 13a failed" << abort(FatalError);

    b2 = cmptMultiply(b1*2.0,pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2.0*b1(i,j,k))
            FatalErrorInFunction << "test 13b failed" << abort(FatalError);

    b2 = cmptMultiply(pTraits<Type>::one,b1);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 13c failed" << abort(FatalError);

    b2 = cmptMultiply(pTraits<Type>::one,b1*2.0);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2.0*b1(i,j,k))
            FatalErrorInFunction << "test 13d failed" << abort(FatalError);


    b2 = cmptDivide(b1,pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != b1(i,j,k))
            FatalErrorInFunction << "test 14a failed" << abort(FatalError);

    b2 = cmptDivide(b1*2.0,pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2.0*b1(i,j,k))
            FatalErrorInFunction << "test 14b failed" << abort(FatalError);

    b2 = cmptDivide(pTraits<Type>::one,(b1+pTraits<Type>::one)());

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(pTraits<Type>::one, b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 14c failed" << abort(FatalError);

    b2 = cmptDivide(pTraits<Type>::one,b1+pTraits<Type>::one);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != cmptDivide(pTraits<Type>::one, b1(i,j,k)+pTraits<Type>::one))
            FatalErrorInFunction << "test 14d failed" << abort(FatalError);

    b2/1.0;
    (1.0*b2)/1.0;
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

    // Let's see if it compiles. Operators are already tested.

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

template<class Type>
void testFaceSpaceFunctions()
{
    const labelVector shape(2,3,4);

    block<FaceSpace<Type>> fs1(shape);
    block<LowerFaceSpace<Type>> ls1(shape);

    block<FaceSpace<scalar>> sfs1(shape);
    block<LowerFaceSpace<scalar>> sls1(shape);

    block<scalar> s1(shape);

    label l = 0;

    forAllBlock(fs1, i, j, k)
    {
        fs1(i,j,k) = pTraits<Type>::one*l++;
        ls1(i,j,k) = pTraits<Type>::one*(l+++1);
        s1(i,j,k) = (l+++1);
        sfs1(i,j,k) = (l+++1);
        sls1(i,j,k) = (l+++1);
    }

    ls1 += fs1;
    ls1 += (1.0*fs1);

    fs1+ls1;
    ls1+fs1;

    (1.0*fs1)+ls1;
    (1.0*ls1)+fs1;

    (1.0*fs1)+(1.0*ls1);
    (1.0*ls1)+(1.0*fs1);

    fs1 += pTraits<Type>::one;
    ls1 += pTraits<Type>::one;

    fs1 += FaceSpace<Type>::one;
    ls1 += FaceSpace<Type>::one;
    ls1 += LowerFaceSpace<Type>::one;

    fs1 + pTraits<Type>::one;
    ls1 + pTraits<Type>::one;

    (1.0*fs1) + pTraits<Type>::one;
    (1.0*ls1) + pTraits<Type>::one;

    fs1 + FaceSpace<Type>::one;
    ls1 + FaceSpace<Type>::one;
    ls1 + LowerFaceSpace<Type>::one;

    (1.0*fs1) + FaceSpace<Type>::one;
    (1.0*ls1) + FaceSpace<Type>::one;
    (1.0*ls1) + LowerFaceSpace<Type>::one;

    pTraits<Type>::one + fs1;
    pTraits<Type>::one + ls1;

    pTraits<Type>::one + (1.0*fs1);
    pTraits<Type>::one + (1.0*ls1);

    FaceSpace<Type>::one + fs1;
    FaceSpace<Type>::one + ls1;
    LowerFaceSpace<Type>::one + ls1;

    FaceSpace<Type>::one + (1.0*fs1);
    FaceSpace<Type>::one + (1.0*ls1);
    LowerFaceSpace<Type>::one + (1.0*ls1);

    ls1 -= fs1;
    ls1 -= (1.0*fs1);

    fs1-ls1;
    ls1-fs1;

    (1.0*fs1)-ls1;
    (1.0*ls1)-fs1;

    fs1-(1.0*ls1);
    ls1-(1.0*fs1);

    fs1 -= pTraits<Type>::one;
    ls1 -= pTraits<Type>::one;

    fs1 -= FaceSpace<Type>::one;
    ls1 -= FaceSpace<Type>::one;
    ls1 -= LowerFaceSpace<Type>::one;

    fs1 - pTraits<Type>::one;
    ls1 - pTraits<Type>::one;

    (1.0*fs1) - pTraits<Type>::one;
    (1.0*ls1) - pTraits<Type>::one;

    fs1 - FaceSpace<Type>::one;
    ls1 - FaceSpace<Type>::one;
    ls1 - LowerFaceSpace<Type>::one;

    (1.0*fs1) - FaceSpace<Type>::one;
    (1.0*ls1) - FaceSpace<Type>::one;
    (1.0*ls1) - LowerFaceSpace<Type>::one;

    pTraits<Type>::one - fs1;
    pTraits<Type>::one - ls1;

    pTraits<Type>::one - (1.0*fs1);
    pTraits<Type>::one - (1.0*ls1);

    FaceSpace<Type>::one - fs1;
    FaceSpace<Type>::one - ls1;
    LowerFaceSpace<Type>::one - ls1;

    FaceSpace<Type>::one - (1.0*fs1);
    FaceSpace<Type>::one - (1.0*ls1);
    LowerFaceSpace<Type>::one - (1.0*ls1);

    fs1 *= (1.0*s1);
    fs1 *= (1.0*sfs1);

    fs1*s1;
    fs1*sfs1;

    (1.0*fs1)*s1;
    (1.0*fs1)*sfs1;

    fs1*(1.0*s1);
    fs1*(1.0*sfs1);

    1.0/s1;
    1.0/sfs1;

    1.0/(1.0*s1);
    1.0/(1.0*sfs1);

    fs1 /= (1.0*s1);
    fs1 /= (1.0*sfs1);

    fs1/s1;
    fs1/sfs1;

    (1.0*fs1)/s1;
    (1.0*fs1)/sfs1;

    fs1/(1.0*s1);
    fs1/(1.0*sfs1);


    ls1 *= (1.0*s1);
    ls1 *= (1.0*sls1);

    ls1*s1;
    ls1*sls1;

    (1.0*ls1)*s1;
    (1.0*ls1)*sls1;

    ls1*(1.0*s1);
    ls1*(1.0*sls1);

    1.0/s1;
    1.0/sls1;

    1.0/(1.0*s1);
    1.0/(1.0*sls1);

    ls1 /= (1.0*s1);
    ls1 /= (1.0*sls1);

    ls1/s1;
    ls1/sls1;

    (1.0*ls1)/s1;
    (1.0*ls1)/sls1;

    ls1/(1.0*s1);
    ls1/(1.0*sls1);
}

void testScalarFunctions()
{
    const labelVector shape(2,3,4);

    block<scalar> b1(shape);

    label l = 0;

    scalar sumsqr = 0;

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = scalar(l+++1);

        sumsqr += Foam::sqr(b1(i,j,k));
    }

    if (sumSqr(b1) != sumsqr)
        FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    if (gSumSqr(b1) != sumsqr)
        FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    scalarBlock b2 = b1/b1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 1)
            FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    b2 = (b1*2.0)/b1;

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 2)
            FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    b2 = b1/(b1*2.0);

    forAllBlock(b2, i, j, k)
        if (b2(i,j,k) != 0.5)
            FatalErrorInFunction << "test 2c failed" << abort(FatalError);
}

void testVectorFunctions()
{
    const labelVector shape(2,3,4);

    block<vector> b1(shape);
    block<vector> b2(shape);

    label l = 0;

    symmTensor b1ss(Zero);

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = vector(l,l+1,l+2); l++;
        b2(i,j,k) = vector(l,l+1,l+2); l++;

        b1ss += Foam::sqr(b1(i,j,k));
    }

    if (sumSqr(b1) != b1ss)
        FatalErrorInFunction << "test 1a failed" << abort(FatalError);

    if (sumSqr(b1*2.0) != 4.0*b1ss)
        FatalErrorInFunction << "test 1b failed" << abort(FatalError);

    if (gSumSqr(b1) != b1ss)
        FatalErrorInFunction << "test 1c failed" << abort(FatalError);

    if (gSumSqr(b1*2.0) != 4.0*b1ss)
        FatalErrorInFunction << "test 1d failed" << abort(FatalError);

    block<tensor> t1 = b1*b2;

    forAllBlock(t1, i, j, k)
        if (t1(i,j,k) != b1(i,j,k)*b2(i,j,k))
            FatalErrorInFunction << "test 2a failed" << abort(FatalError);

    t1 = (b1*2.0)*b2;

    forAllBlock(t1, i, j, k)
        if (t1(i,j,k) != 2.0*b1(i,j,k)*b2(i,j,k))
            FatalErrorInFunction << "test 2b failed" << abort(FatalError);

    t1 = b1*(b2*2.0);

    forAllBlock(t1, i, j, k)
        if (t1(i,j,k) != 2.0*b1(i,j,k)*b2(i,j,k))
            FatalErrorInFunction << "test 2c failed" << abort(FatalError);

    t1 = (b1*2.0)*(b2*2.0);

    forAllBlock(t1, i, j, k)
        if (t1(i,j,k) != 4.0*b1(i,j,k)*b2(i,j,k))
            FatalErrorInFunction << "test 2d failed" << abort(FatalError);


    scalarBlock s1 = b1 & b2;

    forAllBlock(s1, i, j, k)
        if (s1(i,j,k) != (b1(i,j,k) & b2(i,j,k)))
            FatalErrorInFunction << "test 3a failed" << abort(FatalError);

    s1 = (b1*2.0) & b2;

    forAllBlock(s1, i, j, k)
        if (s1(i,j,k) != 2.0*(b1(i,j,k) & b2(i,j,k)))
            FatalErrorInFunction << "test 3b failed" << abort(FatalError);

    s1 = b1 & (b2*2.0);

    forAllBlock(s1, i, j, k)
        if (s1(i,j,k) != 2.0*(b1(i,j,k) & b2(i,j,k)))
            FatalErrorInFunction << "test 3c failed" << abort(FatalError);

    s1 = (b1*2.0) & (b2*2.0);

    forAllBlock(s1, i, j, k)
        if (s1(i,j,k) != 4.0*(b1(i,j,k) & b2(i,j,k)))
            FatalErrorInFunction << "test 3d failed" << abort(FatalError);


    vectorBlock v1 = b1 ^ b2;

    forAllBlock(s1, i, j, k)
        if (v1(i,j,k) != (b1(i,j,k) ^ b2(i,j,k)))
            FatalErrorInFunction << "test 4a failed" << abort(FatalError);

    v1 = (b1*2.0) ^ b2;

    forAllBlock(s1, i, j, k)
        if (v1(i,j,k) != 2.0*(b1(i,j,k) ^ b2(i,j,k)))
            FatalErrorInFunction << "test 4b failed" << abort(FatalError);

    v1 = b1 ^ (b2*2.0);

    forAllBlock(s1, i, j, k)
        if (v1(i,j,k) != 2.0*(b1(i,j,k) ^ b2(i,j,k)))
            FatalErrorInFunction << "test 4c failed" << abort(FatalError);

    v1 = (b1*2.0) ^ (b2*2.0);

    forAllBlock(s1, i, j, k)
        if (v1(i,j,k) != 4.0*(b1(i,j,k) ^ b2(i,j,k)))
            FatalErrorInFunction << "test 4d failed" << abort(FatalError);
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

    // Macro function to check if collective operator results in the same as the
    // operator applied to individual entries

    #define TEST(RET, B1, B2, OP, NAME)                                         \
    {                                                                           \
        block<RET> r = B1 OP B2;                                                \
                                                                                \
        forAllBlock(r, i, j, k)                                                 \
            if (r(i,j,k) != (B1(i,j,k) OP B2(i,j,k)))                           \
                FatalErrorInFunction                                            \
                    << "test " NAME "a failed"                                  \
                    << abort(FatalError);                                       \
                                                                                \
        r = (B1*2.0) OP B2;                                                     \
                                                                                \
        forAllBlock(r, i, j, k)                                                 \
            if (r(i,j,k) != (2.0*B1(i,j,k) OP B2(i,j,k)))                       \
                FatalErrorInFunction                                            \
                    << "test " NAME "b failed"                                  \
                    << abort(FatalError);                                       \
                                                                                \
        r = B1 OP (B2*2.0);                                                     \
                                                                                \
        forAllBlock(r, i, j, k)                                                 \
            if (r(i,j,k) != 2.0*(B1(i,j,k) OP B2(i,j,k)))                       \
                FatalErrorInFunction                                            \
                    << "test " NAME "c failed"                                  \
                    << abort(FatalError);                                       \
                                                                                \
        r = (B1*2.0) OP (B2*2.0);                                               \
                                                                                \
        forAllBlock(r, i, j, k)                                                 \
            if (r(i,j,k) != 4.0*(B1(i,j,k) OP B2(i,j,k)))                       \
                FatalErrorInFunction                                            \
                    << "test " NAME "d failed"                                  \
                    << abort(FatalError);                                       \
    }

    // Commented lines don't work because of OpenFOAM bugs or lack of
    // implementation

    TEST(tensor, b1, b1, &, "1")
    TEST(tensor, b1, b2, &, "2")
    TEST(tensor, b1, b3, &, "3")
    TEST(tensor, b1, b4, &, "4")
    TEST(tensor, b2, b1, &, "5")
    TEST(tensor, b3, b1, &, "6")
    TEST(tensor, b4, b1, &, "7")
    TEST(tensor, b2, b2, &, "8")
    TEST(symmTensor, b2, b3, &, "9")
    //TEST(tensor, b2, b4, &, "10")
    TEST(symmTensor, b3, b2, &, "11")
    // TEST(tensor, b4, b2, &, "12")
    TEST(sphericalTensor, b3, b3, &, "13")
    // TEST(tensor, b3, b4, &, "14")
    // TEST(tensor, b4, b3, &, "15")
    // TEST(tensor, b4, b4, &, "16")

    TEST(scalar, b1, b1, &&, "17")
    TEST(scalar, b1, b2, &&, "18")
    TEST(scalar, b1, b3, &&, "19")
    // TEST(scalar, b1, b4, &&, "20")
    TEST(scalar, b2, b1, &&, "21")
    TEST(scalar, b3, b1, &&, "22")
    // TEST(scalar, b4, b1, &&, "23")
    TEST(scalar, b2, b2, &&, "24")
    TEST(scalar, b2, b3, &&, "25")
    // TEST(scalar, b2, b4, &&, "26")
    TEST(scalar, b3, b2, &&, "27")
    // TEST(scalar, b4, b2, &&, "28")
    TEST(scalar, b3, b3, &&, "29")
    // TEST(scalar, b3, b4, &&, "30")
    // TEST(scalar, b4, b3, &&, "31")
    TEST(scalar, b4, b4, &&, "32")

    TEST(vector, b1, v1, &, "33")
    TEST(vector, b2, v1, &, "34")
    TEST(vector, b3, v1, &, "35")
    TEST(vector, b4, v1, &, "36")
    TEST(vector, v1, b1, &, "37")
    TEST(vector, v1, b2, &, "38")
    TEST(vector, v1, b3, &, "39")
    TEST(vector, v1, b4, &, "40")

    #undef TEST
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    testConstructors<label>();
    testConstructors<scalar>();

    testConstructors<vector>();
    testConstructors<tensor>();
    testConstructors<symmTensor>();
    testConstructors<sphericalTensor>();
    testConstructors<diagTensor>();

    testConstructors<faceScalar>();
    testConstructors<lowerFaceScalar>();
    testConstructors<edgeScalar>();
    testConstructors<vertexScalar>();

    testConstructors<faceVector>();
    testConstructors<lowerFaceVector>();
    testConstructors<edgeVector>();
    testConstructors<vertexVector>();

    testConstructors<stencil>();
    testConstructors<diagStencil>();


    testIndexing<label>();
    testIndexing<scalar>();

    testIndexing<vector>();
    testIndexing<tensor>();
    testIndexing<symmTensor>();
    testIndexing<sphericalTensor>();
    testIndexing<diagTensor>();

    testIndexing<faceScalar>();
    testIndexing<lowerFaceScalar>();
    testIndexing<edgeScalar>();
    testIndexing<vertexScalar>();

    testIndexing<faceVector>();
    testIndexing<lowerFaceVector>();
    testIndexing<edgeVector>();
    testIndexing<vertexVector>();

    testIndexing<stencil>();
    testIndexing<diagStencil>();


    testTransformations<label>();
    testTransformations<scalar>();

    testTransformations<vector>();
    testTransformations<tensor>();
    testTransformations<symmTensor>();
    testTransformations<sphericalTensor>();
    testTransformations<diagTensor>();

    testTransformations<faceScalar>();
    testTransformations<lowerFaceScalar>();
    testTransformations<edgeScalar>();
    testTransformations<vertexScalar>();

    testTransformations<faceVector>();
    testTransformations<lowerFaceVector>();
    testTransformations<edgeVector>();
    testTransformations<vertexVector>();

    testTransformations<stencil>();
    testTransformations<diagStencil>();


    testMemberOperators<label>();
    testMemberOperators<scalar>();

    testMemberOperators<vector>();
    testMemberOperators<tensor>();
    testMemberOperators<symmTensor>();
    testMemberOperators<sphericalTensor>();
    testMemberOperators<diagTensor>();

    testMemberOperators<faceScalar>();
    testMemberOperators<lowerFaceScalar>();
    testMemberOperators<edgeScalar>();
    testMemberOperators<vertexScalar>();

    testMemberOperators<faceVector>();
    testMemberOperators<lowerFaceVector>();
    testMemberOperators<edgeVector>();
    testMemberOperators<vertexVector>();

    testMemberOperators<stencil>();
    testMemberOperators<diagStencil>();


    testPrimitiveFunctions<scalar>();

    testPrimitiveFunctions<vector>();
    testPrimitiveFunctions<tensor>();
    testPrimitiveFunctions<symmTensor>();
    testPrimitiveFunctions<sphericalTensor>();
    testPrimitiveFunctions<diagTensor>();

    testPrimitiveFunctions<faceScalar>();
    testPrimitiveFunctions<lowerFaceScalar>();
    testPrimitiveFunctions<edgeScalar>();
    testPrimitiveFunctions<vertexScalar>();

    testPrimitiveFunctions<faceVector>();
    testPrimitiveFunctions<lowerFaceVector>();
    testPrimitiveFunctions<edgeVector>();
    testPrimitiveFunctions<vertexVector>();


    testVectorSpaceFunctions<vector>();
    testVectorSpaceFunctions<tensor>();
    testVectorSpaceFunctions<symmTensor>();
    testVectorSpaceFunctions<sphericalTensor>();
    testVectorSpaceFunctions<diagTensor>();


    testStencilFunctions<stencil>();
    testStencilFunctions<diagStencil>();


    testFaceSpaceFunctions<scalar>();
    testFaceSpaceFunctions<vector>();

    testScalarFunctions();
    testVectorFunctions();
    testTensorFunctions();
}
