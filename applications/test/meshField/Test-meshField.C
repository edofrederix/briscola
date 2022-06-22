#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testConstructors(const fvMesh& fvMsh)
{
    // Field from name and mesh

    meshField<Type,MeshType> m1("m1", fvMsh);

    // Read

    meshField<Type,MeshType> m2
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ
    );

    // Copy from m2 with same name

    meshField<Type,MeshType> m3a(m2);
    meshField<Type,MeshType> m3b(m2, true);
    meshField<Type,MeshType> m3c(m2, true, true);

    // Copy from m2 with new name

    meshField<Type,MeshType> m4a("m4a", m2);
    meshField<Type,MeshType> m4b("m4b", m2, true);
    meshField<Type,MeshType> m4c("m4c", m2, true, true);

    // Copy from tmp of m2 with same name

    meshField<Type,MeshType> m5a(2.0*m2);
    meshField<Type,MeshType> m5b(2.0*m2, true);
    meshField<Type,MeshType> m5c(2.0*m2, true, true);

    // Copy from tmp m2 with new name

    meshField<Type,MeshType> m6a("m6a", 2.0*m2);
    meshField<Type,MeshType> m6b("m6b", 2.0*m2, true);
    meshField<Type,MeshType> m6c("m6c", 2.0*m2, true, true);
}

template<class Type, class MeshType>
void testIndexing(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ
    );

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
    }

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        if (m1[l][d](i,j,k) != m1[l][d](labelVector(i,j,k)))
            FatalErrorInFunction << "test 2 failed" << abort(FatalError);
    }
}

template<class Type, class MeshType>
void testMemberOperators(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1[l][d](i,j,k) = c++;
    }

    m1 = Zero;
    m1 = pTraits<Type>::one*2.0;
    m1 = m2;
    m1 = (2.0*m2);

    m1 += m2;
    m1 += (2.0*m2);

    m1 -= m2;
    m1 -= (2.0*m2);

    m1 *= s1;
    m1 *= (2.0*s1);

    m1 /= s1;
    m1 /= (2.0*s1);

    m1 += pTraits<Type>::one;
    m1 -= pTraits<Type>::one;

    m1 *= scalar(2.0);
    m1 /= scalar(2.0);
}

template<class Type, class MeshType>
void testPrimitiveFunctions(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1[l][d](i,j,k) = c++;
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
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
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
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1[l][d](i,j,k) = c++;
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
    meshField<scalar,MeshType> m1("m1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = scalar(c+++1);
    }

    m1/m1;
    (m1*2.0)/m1;
    m1/(m1*2.0);
}

template<class MeshType>
void testVectorFunctions(const fvMesh& fvMsh)
{
    meshField<vector,MeshType> m1("m1", fvMsh);
    meshField<vector,MeshType> m2("m2", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<vector>::one*c++;
        m2[l][d](i,j,k) = pTraits<vector>::one*c++;
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
    meshField<tensor,MeshType> m1("m1", fvMsh);
    meshField<symmTensor,MeshType> m2("m2", fvMsh);
    meshField<sphericalTensor,MeshType> b3("b3", fvMsh);
    meshField<diagTensor,MeshType> b4("b4", fvMsh);
    meshField<vector,MeshType> v1("v1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<tensor>::one*c++;
        m2[l][d](i,j,k) = pTraits<symmTensor>::one*c++;
        b3[l][d](i,j,k) = pTraits<sphericalTensor>::one*c++;
        b4[l][d](i,j,k) = pTraits<diagTensor>::one*c++;
        v1[l][d](i,j,k) = vector(l,l+1,l+2); c++;
    }

    m1 & m1;
    (m1*2.0) & m1;
    m1 & (m1*2.0);
    (m1*2.0) & (m1*2.0);

    m1 & m2;
    (m1*2.0) & m2;
    m1 & (m2*2.0);
    (m1*2.0) & (m2*2.0);

    m1 && m1;
    (m1*2.0) && m1;
    m1 && (m1*2.0);
    (m1*2.0) && (m1*2.0);

    m1 && m2;
    (m1*2.0) && m2;
    m1 && (m2*2.0);
    (m1*2.0) && (m2*2.0);

    m1 & v1;
    (m1*2.0) & v1;
    v1 & (m1*2.0);
    (v1*2.0) & (m1*2.0);

    m2 & v1;
    (m2*2.0) & v1;
    v1 & (m2*2.0);
    (v1*2.0) & (m2*2.0);

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
    testConstructors<hexScalar,colocated>(fvMsh);
    testConstructors<vector,colocated>(fvMsh);
    testConstructors<hexVector,colocated>(fvMsh);
    testConstructors<tensor,colocated>(fvMsh);
    testConstructors<symmTensor,colocated>(fvMsh);
    testConstructors<sphericalTensor,colocated>(fvMsh);
    testConstructors<diagTensor,colocated>(fvMsh);
    testConstructors<stencil,colocated>(fvMsh);
    testConstructors<symmStencil,colocated>(fvMsh);
    testConstructors<diagStencil,colocated>(fvMsh);

    testIndexing<label,colocated>(fvMsh);
    testIndexing<scalar,colocated>(fvMsh);
    testIndexing<hexScalar,colocated>(fvMsh);
    testIndexing<vector,colocated>(fvMsh);
    testIndexing<hexVector,colocated>(fvMsh);
    testIndexing<tensor,colocated>(fvMsh);
    testIndexing<symmTensor,colocated>(fvMsh);
    testIndexing<sphericalTensor,colocated>(fvMsh);
    testIndexing<diagTensor,colocated>(fvMsh);
    testIndexing<stencil,colocated>(fvMsh);
    testIndexing<symmStencil,colocated>(fvMsh);
    testIndexing<diagStencil,colocated>(fvMsh);

    testMemberOperators<label,colocated>(fvMsh);
    testMemberOperators<scalar,colocated>(fvMsh);
    testMemberOperators<hexScalar,colocated>(fvMsh);
    testMemberOperators<vector,colocated>(fvMsh);
    testMemberOperators<hexVector,colocated>(fvMsh);
    testMemberOperators<tensor,colocated>(fvMsh);
    testMemberOperators<symmTensor,colocated>(fvMsh);
    testMemberOperators<sphericalTensor,colocated>(fvMsh);
    testMemberOperators<diagTensor,colocated>(fvMsh);
    testMemberOperators<stencil,colocated>(fvMsh);
    testMemberOperators<symmStencil,colocated>(fvMsh);
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
    testStencilFunctions<symmStencil,colocated>(fvMsh);
    testStencilFunctions<diagStencil,colocated>(fvMsh);

    testScalarFunctions<colocated>(fvMsh);
    testVectorFunctions<colocated>(fvMsh);
    testTensorFunctions<colocated>(fvMsh);

    // Staggered

    testConstructors<label,staggered>(fvMsh);
    testConstructors<scalar,staggered>(fvMsh);
    testConstructors<hexScalar,staggered>(fvMsh);
    testConstructors<vector,staggered>(fvMsh);
    testConstructors<hexVector,staggered>(fvMsh);
    testConstructors<tensor,staggered>(fvMsh);
    testConstructors<symmTensor,staggered>(fvMsh);
    testConstructors<sphericalTensor,staggered>(fvMsh);
    testConstructors<diagTensor,staggered>(fvMsh);
    testConstructors<stencil,staggered>(fvMsh);
    testConstructors<symmStencil,staggered>(fvMsh);
    testConstructors<diagStencil,staggered>(fvMsh);

    testIndexing<label,staggered>(fvMsh);
    testIndexing<scalar,staggered>(fvMsh);
    testIndexing<hexScalar,staggered>(fvMsh);
    testIndexing<vector,staggered>(fvMsh);
    testIndexing<hexVector,staggered>(fvMsh);
    testIndexing<tensor,staggered>(fvMsh);
    testIndexing<symmTensor,staggered>(fvMsh);
    testIndexing<sphericalTensor,staggered>(fvMsh);
    testIndexing<diagTensor,staggered>(fvMsh);
    testIndexing<stencil,staggered>(fvMsh);
    testIndexing<symmStencil,staggered>(fvMsh);
    testIndexing<diagStencil,staggered>(fvMsh);

    testMemberOperators<label,staggered>(fvMsh);
    testMemberOperators<scalar,staggered>(fvMsh);
    testMemberOperators<hexScalar,staggered>(fvMsh);
    testMemberOperators<vector,staggered>(fvMsh);
    testMemberOperators<hexVector,staggered>(fvMsh);
    testMemberOperators<tensor,staggered>(fvMsh);
    testMemberOperators<symmTensor,staggered>(fvMsh);
    testMemberOperators<sphericalTensor,staggered>(fvMsh);
    testMemberOperators<diagTensor,staggered>(fvMsh);
    testMemberOperators<stencil,staggered>(fvMsh);
    testMemberOperators<symmStencil,staggered>(fvMsh);
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
    testStencilFunctions<symmStencil,staggered>(fvMsh);
    testStencilFunctions<diagStencil,staggered>(fvMsh);

    testScalarFunctions<staggered>(fvMsh);
    testVectorFunctions<staggered>(fvMsh);
    testTensorFunctions<staggered>(fvMsh);
}
