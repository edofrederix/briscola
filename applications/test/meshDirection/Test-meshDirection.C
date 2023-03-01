#include "arguments.H"
#include "Time.H"

#include "colocatedDirections.H"
#include "staggeredDirections.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testConstructors(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0);
    meshDirection<Type,MeshType> m2(fvMsh,0,0,Zero);
    meshDirection<Type,MeshType> m3(fvMsh,0,0,pTraits<Type>::one);

    meshDirection<Type,MeshType> m4(m3);
    meshDirection<Type,MeshType> m5(m3, Zero);
    meshDirection<Type,MeshType> m6(m3, pTraits<Type>::one);

    meshDirection<Type,MeshType> m7(m3*2.0);
    meshDirection<Type,MeshType> m8(m3*2.0, Zero);
    meshDirection<Type,MeshType> m9(m3*2.0, pTraits<Type>::one);
}

template<class Type, class MeshType>
void testIndexing(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0,pTraits<Type>::one);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<Type>::one*l++;
    }

    forAllCells(m1, i, j, k)
    {
        if (m1(i,j,k) != m1(labelVector(i,j,k)))
            FatalErrorInFunction << "test 2 failed" << abort(FatalError);
    }
}

template<class Type, class MeshType>
void testMemberOperators(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0);
    meshDirection<Type,MeshType> m2(fvMsh,0,0);

    meshDirection<scalar,MeshType> s1(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<Type>::one*l++;
        m2(i,j,k) = pTraits<Type>::one*l++;

        s1(i,j,k) = l++;
    }

    m1 = Zero;
    m1 = pTraits<Type>::one*2.0;
    m1 = m2;
    m1 = m2*2;

    m1 += m2;
    m1 += m2*2;

    m1 -= m2;
    m1 -= m2*2;

    m1 *= s1;
    m1 *= s1*2;

    m1 /= s1;
    m1 /= s1*2;

    m1 += pTraits<Type>::one;
    m1 -= pTraits<Type>::one;

    m1 + pTraits<Type>::one;
    pTraits<Type>::one + m1;

    m1 - pTraits<Type>::one;
    pTraits<Type>::one - m1;

    m1 *= scalar(2.0);
    m1 /= scalar(2.0);
}

template<class Type, class MeshType>
void testPrimitiveFunctions(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0);
    meshDirection<Type,MeshType> m2(fvMsh,0,0);

    meshDirection<scalar,MeshType> s1(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<Type>::one*l++;
        m2(i,j,k) = pTraits<Type>::one*l++;

        s1(i,j,k) = l++;
    }

    mag(m1);
    mag(m1*2.0);

    max(m1);
    max(m1*2.0);
    gMax(m1);
    gMax(m1*2.0);

    min(m1);
    min(m1*2.0);
    gMin(m1);
    gMin(m1*2.0);

    sum(m1);
    sum(m1*2.0);
    gSum(m1);
    gSum(m1*2.0);

    average(m1);
    average(m1*2.0);
    gAverage(m1);
    gAverage(m1*2.0);

    sumProd(m1, m1);
    sumProd(m1*2.0, m1);
    sumProd(m1, m1*2.0);
    sumProd(m1*2.0, m1*2.0);

    max(m1,m2);
    max(m1*2.0,m2);
    max(m1,m2*2.0);
    max(m1*2.0,m2*2.0);

    min(m1,m2);
    min(m1*2.0,m2);
    min(m1,m2*2.0);
    min(m1*2.0,m2*2.0);

    max(m1,pTraits<Type>::one);
    max(m1*2.0,pTraits<Type>::one);
    max(pTraits<Type>::one,m2);
    max(pTraits<Type>::one,m2*2.0);

    min(m1,pTraits<Type>::one);
    min(m1*2.0,pTraits<Type>::one);
    min(pTraits<Type>::one,m2);
    min(pTraits<Type>::one,m2*2.0);

    -m2;

    m1*s1;
    (m1*2.0)*s1;
    s1*m1;
    s1*(m1*2.0);

    m1/s1;
    (m1*2.0)/s1;

    m1+m2;

    m1-m2;
}

template<class Type, class MeshType>
void testVectorSpaceFunctions(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0);
    meshDirection<Type,MeshType> m2(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<Type>::one*l++;
        m2(i,j,k) = pTraits<Type>::one*l++;
    }

    cmptMax(m1);
    cmptMax(m1*2.0);

    cmptMin(m1);
    cmptMin(m1*2.0);

    cmptAv(m1);
    cmptAv(m1*2.0);

    cmptMag(m1);
    cmptMag(m1*2.0);

    maxMagSqr(m1);
    maxMagSqr(m1*2.0);
    gMaxMagSqr(m1);
    gMaxMagSqr(m1*2.0);

    minMagSqr(m1);
    minMagSqr(m1*2.0);
    gMinMagSqr(m1);
    gMinMagSqr(m1*2.0);

    sumMag(m1);
    sumMag(m1*2.0);
    gSumMag(m1);
    gSumMag(m1*2.0);

    sumCmptProd(m1, m1);
    sumCmptProd(m1*2.0, m1);
    sumCmptProd(m1, m1*2.0);
    sumCmptProd(m1*2.0, m1*2.0);

    sumCmptMag(m1);
    sumCmptMag(m1*2.0);
    gSumCmptMag(m1);
    gSumCmptMag(m1*2.0);

    cmptMultiply(m1,m2);
    cmptMultiply(m1*2.0,m2);
    cmptMultiply(m1,m2*2.0);
    cmptMultiply(m1*2.0,m2*2.0);

    cmptDivide(m1,m2);
    cmptDivide(m1*2.0,m2);
    cmptDivide(m1,m2*2.0);
    cmptDivide(m1*2.0,m2*2.0);

    cmptMultiply(m1,pTraits<Type>::one);
    cmptMultiply(m1*2.0,pTraits<Type>::one);
    cmptMultiply(pTraits<Type>::one,m2);
    cmptMultiply(pTraits<Type>::one,m2*2.0);

    cmptDivide(m1,pTraits<Type>::one);
    cmptDivide(m1*2.0,pTraits<Type>::one);
    cmptDivide(pTraits<Type>::one,m2);
    cmptDivide(pTraits<Type>::one,m2*2.0);
}

template<class Type, class MeshType>
void testStencilFunctions(const fvMesh& fvMsh)
{
    meshDirection<Type,MeshType> m1(fvMsh,0,0);
    meshDirection<Type,MeshType> m2(fvMsh,0,0);

    meshDirection<scalar,MeshType> s1(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<Type>::one*l++;
        m2(i,j,k) = pTraits<Type>::one*l++;

        s1(i,j,k) = l++;
    }

    -m1;

    m1+m2;
    (2.0*m1)+m2;
    m1+(2.0*m2);
    (2.0*m1)+(2.0*m2);

    m1-m2;
    (2.0*m1)-m2;
    m1-(2.0*m2);
    (2.0*m1)-(2.0*m2);

    m1-m2;
    (2.0*m1)-m2;
    m1-(2.0*m2);
    (2.0*m1)-(2.0*m2);

    m1*s1;
    s1*m1;
    (2.0*m1)*s1;
    s1*(2.0*m1);

    m1*(2.0*s1);
    (2.0*s1)*m1;
    (2.0*m1)*(2.0*s1);
    (2.0*s1)*(2.0*m1);

    m1/s1;
    (2.0*m1)/s1;

    m1/(2.0*s1);
    (2.0*m1)/(2.0*s1);
}

template<class MeshType>
void testScalarFunctions(const fvMesh& fvMsh)
{
    meshDirection<scalar,MeshType> m1(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = scalar(l+++1);
    }

    m1/m1;
    (m1*2.0)/m1;
    m1/(m1*2.0);
}

template<class MeshType>
void testVectorFunctions(const fvMesh& fvMsh)
{
    meshDirection<vector,MeshType> m1(fvMsh,0,0);
    meshDirection<vector,MeshType> m2(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = vector(l,l+1,l+2); l++;
        m2(i,j,k) = vector(l,l+1,l+2); l++;
    }

    m1*m2;
    (m1*2.0)*m2;
    m1*(m2*2.0);
    (m1*2.0)*(m2*2.0);

    m1 & m2;
    (m1*2.0) & m2;
    m1 & (m2*2.0);
    (m1*2.0) & (m2*2.0);

    m1 ^ m2;
    (m1*2.0) ^ m2;
    m1 ^ (m2*2.0);
    (m1*2.0) ^ (m2*2.0);
}

template<class MeshType>
void testTensorFunctions(const fvMesh& fvMsh)
{
    meshDirection<tensor,MeshType> m1(fvMsh,0,0);
    meshDirection<symmTensor,MeshType> m2(fvMsh,0,0);
    meshDirection<sphericalTensor,MeshType> m3(fvMsh,0,0);
    meshDirection<diagTensor,MeshType> m4(fvMsh,0,0);
    meshDirection<vector,MeshType> v1(fvMsh,0,0);

    label l = 0;

    forAllCells(m1, i, j, k)
    {
        m1(i,j,k) = pTraits<tensor>::one*l++;
        m2(i,j,k) = pTraits<symmTensor>::one*l++;
        m3(i,j,k) = pTraits<sphericalTensor>::one*l++;
        m4(i,j,k) = pTraits<diagTensor>::one*l++;
        v1(i,j,k) = vector(l,l+1,l+2); l++;
    }

    // Commented combinations could exist, but are not defined as primitives in
    // OpenFOAM

    m1 & m1;
    (m1*2.0) & m1;
    m1 & (m1*2.0);
    (m1*2.0) & (m1*2.0);

    m1 & m2;
    (m1*2.0) & m2;
    m1 & (m2*2.0);
    (m1*2.0) & (m2*2.0);

    m1 & m3;
    (m1*2.0) & m3;
    m1 & (m3*2.0);
    (m1*2.0) & (m3*2.0);

    m1 & m4;
    (m1*2.0) & m4;
    m1 & (m4*2.0);
    (m1*2.0) & (m4*2.0);

    m1 & v1;
    (m1*2.0) & v1;
    m1 & (v1*2.0);
    (m1*2.0) & (v1*2.0);

    m2 & m1;
    (m2*2.0) & m1;
    m2 & (m1*2.0);
    (m2*2.0) & (m1*2.0);

    m3 & m1;
    (m3*2.0) & m1;
    m3 & (m1*2.0);
    (m3*2.0) & (m1*2.0);

    m4 & m1;
    (m4*2.0) & m1;
    m4 & (m1*2.0);
    (m4*2.0) & (m1*2.0);

    v1 & m1;
    (v1*2.0) & m1;
    v1 & (m1*2.0);
    (v1*2.0) & (m1*2.0);

    m1 && m1;
    (m1*2.0) && m1;
    m1 && (m1*2.0);
    (m1*2.0) && (m1*2.0);

    m1 && m2;
    (m1*2.0) && m2;
    m1 && (m2*2.0);
    (m1*2.0) && (m2*2.0);

    m1 && m3;
    (m1*2.0) && m3;
    m1 && (m3*2.0);
    (m1*2.0) && (m3*2.0);

    // m1 && m4;
    // (m1*2.0) && m4;
    // m1 && (m4*2.0);
    // (m1*2.0) && (m4*2.0);

    // m1 && v1;
    // (m1*2.0) && v1;
    // m1 && (v1*2.0);
    // (m1*2.0) && (v1*2.0);

    m2 && m1;
    (m2*2.0) && m1;
    m2 && (m1*2.0);
    (m2*2.0) && (m1*2.0);

    m3 && m1;
    (m3*2.0) && m1;
    m3 && (m1*2.0);
    (m3*2.0) && (m1*2.0);

    // m4 && m1;
    // (m4*2.0) && m1;
    // m4 && (m1*2.0);
    // (m4*2.0) && (m1*2.0);

    // v1 && m1;
    // (v1*2.0) && m1;
    // v1 && (m1*2.0);
    // (v1*2.0) && (m1*2.0);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    // Colocated

    testConstructors<label,colocated>(fvMsh);
    testConstructors<scalar,colocated>(fvMsh);

    testConstructors<vector,colocated>(fvMsh);
    testConstructors<tensor,colocated>(fvMsh);
    testConstructors<symmTensor,colocated>(fvMsh);
    testConstructors<sphericalTensor,colocated>(fvMsh);
    testConstructors<diagTensor,colocated>(fvMsh);

    testConstructors<faceScalar,colocated>(fvMsh);
    testConstructors<faceVector,colocated>(fvMsh);

    testConstructors<stencil,colocated>(fvMsh);
    testConstructors<diagStencil,colocated>(fvMsh);


    testIndexing<label,colocated>(fvMsh);
    testIndexing<scalar,colocated>(fvMsh);

    testIndexing<vector,colocated>(fvMsh);
    testIndexing<tensor,colocated>(fvMsh);
    testIndexing<symmTensor,colocated>(fvMsh);
    testIndexing<sphericalTensor,colocated>(fvMsh);
    testIndexing<diagTensor,colocated>(fvMsh);

    testIndexing<faceScalar,colocated>(fvMsh);
    testIndexing<faceVector,colocated>(fvMsh);

    testIndexing<stencil,colocated>(fvMsh);
    testIndexing<diagStencil,colocated>(fvMsh);


    testMemberOperators<label,colocated>(fvMsh);
    testMemberOperators<scalar,colocated>(fvMsh);

    testMemberOperators<vector,colocated>(fvMsh);
    testMemberOperators<tensor,colocated>(fvMsh);
    testMemberOperators<symmTensor,colocated>(fvMsh);
    testMemberOperators<sphericalTensor,colocated>(fvMsh);
    testMemberOperators<diagTensor,colocated>(fvMsh);

    testMemberOperators<faceScalar,colocated>(fvMsh);
    testMemberOperators<faceVector,colocated>(fvMsh);

    testMemberOperators<stencil,colocated>(fvMsh);
    testMemberOperators<diagStencil,colocated>(fvMsh);


    testPrimitiveFunctions<label,colocated>(fvMsh);
    testPrimitiveFunctions<scalar,colocated>(fvMsh);

    testPrimitiveFunctions<vector,colocated>(fvMsh);
    testPrimitiveFunctions<tensor,colocated>(fvMsh);
    testPrimitiveFunctions<symmTensor,colocated>(fvMsh);
    testPrimitiveFunctions<sphericalTensor,colocated>(fvMsh);
    testPrimitiveFunctions<diagTensor,colocated>(fvMsh);

    testVectorSpaceFunctions<vector,colocated>(fvMsh);
    testVectorSpaceFunctions<tensor,colocated>(fvMsh);
    testVectorSpaceFunctions<symmTensor,colocated>(fvMsh);
    testVectorSpaceFunctions<sphericalTensor,colocated>(fvMsh);
    testVectorSpaceFunctions<diagTensor,colocated>(fvMsh);

    testStencilFunctions<stencil,colocated>(fvMsh);
    testStencilFunctions<diagStencil,colocated>(fvMsh);

    testScalarFunctions<colocated>(fvMsh);
    testVectorFunctions<colocated>(fvMsh);
    testTensorFunctions<colocated>(fvMsh);

    // Staggered

    if (fvMsh.structured())
    {
        testConstructors<label,staggered>(fvMsh);
        testConstructors<scalar,staggered>(fvMsh);
        testConstructors<vector,staggered>(fvMsh);

        testConstructors<tensor,staggered>(fvMsh);
        testConstructors<symmTensor,staggered>(fvMsh);
        testConstructors<sphericalTensor,staggered>(fvMsh);
        testConstructors<diagTensor,staggered>(fvMsh);

        testConstructors<faceScalar,staggered>(fvMsh);
        testConstructors<faceVector,staggered>(fvMsh);

        testConstructors<stencil,staggered>(fvMsh);
        testConstructors<diagStencil,staggered>(fvMsh);


        testIndexing<label,staggered>(fvMsh);
        testIndexing<scalar,staggered>(fvMsh);

        testIndexing<vector,staggered>(fvMsh);
        testIndexing<tensor,staggered>(fvMsh);
        testIndexing<symmTensor,staggered>(fvMsh);
        testIndexing<sphericalTensor,staggered>(fvMsh);
        testIndexing<diagTensor,staggered>(fvMsh);

        testIndexing<faceScalar,staggered>(fvMsh);
        testIndexing<faceVector,staggered>(fvMsh);

        testIndexing<stencil,staggered>(fvMsh);
        testIndexing<diagStencil,staggered>(fvMsh);


        testMemberOperators<label,staggered>(fvMsh);
        testMemberOperators<scalar,staggered>(fvMsh);

        testMemberOperators<vector,staggered>(fvMsh);
        testMemberOperators<tensor,staggered>(fvMsh);
        testMemberOperators<symmTensor,staggered>(fvMsh);
        testMemberOperators<sphericalTensor,staggered>(fvMsh);
        testMemberOperators<diagTensor,staggered>(fvMsh);

        testMemberOperators<faceScalar,staggered>(fvMsh);
        testMemberOperators<faceVector,staggered>(fvMsh);

        testMemberOperators<stencil,staggered>(fvMsh);
        testMemberOperators<diagStencil,staggered>(fvMsh);


        testPrimitiveFunctions<label,staggered>(fvMsh);
        testPrimitiveFunctions<scalar,staggered>(fvMsh);

        testPrimitiveFunctions<vector,staggered>(fvMsh);
        testPrimitiveFunctions<tensor,staggered>(fvMsh);
        testPrimitiveFunctions<symmTensor,staggered>(fvMsh);
        testPrimitiveFunctions<sphericalTensor,staggered>(fvMsh);
        testPrimitiveFunctions<diagTensor,staggered>(fvMsh);


        testVectorSpaceFunctions<vector,staggered>(fvMsh);
        testVectorSpaceFunctions<tensor,staggered>(fvMsh);
        testVectorSpaceFunctions<symmTensor,staggered>(fvMsh);
        testVectorSpaceFunctions<sphericalTensor,staggered>(fvMsh);
        testVectorSpaceFunctions<diagTensor,staggered>(fvMsh);


        testStencilFunctions<stencil,staggered>(fvMsh);
        testStencilFunctions<diagStencil,staggered>(fvMsh);


        testScalarFunctions<staggered>(fvMsh);
        testVectorFunctions<staggered>(fvMsh);
        testTensorFunctions<staggered>(fvMsh);
    }
}
