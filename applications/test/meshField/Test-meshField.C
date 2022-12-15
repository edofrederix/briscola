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

    meshField<Type,MeshType> m5a(2*m2);
    meshField<Type,MeshType> m5b(2*m2, true);
    meshField<Type,MeshType> m5c(2*m2, true, true);

    // Copy from tmp m2 with new name

    meshField<Type,MeshType> m6a("m6a", 2*m2);
    meshField<Type,MeshType> m6b("m6b", 2*m2, true);
    meshField<Type,MeshType> m6c("m6c", 2*m2, true, true);
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
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
    }

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        if (m1[l][d](i,j,k) != m1[l][d](labelVector(i,j,k)))
            FatalErrorInFunction << "test 1 failed" << abort(FatalError);
    }
}

template<class Type, class MeshType>
void testMemberOperators(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);
        s1[l][d](i,j,k) = scalar(l+d+i+j+k+1);
    }

    m2 = m1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1a failed" << abort(FatalError);

    m2 = 1.0*m1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1b failed" << abort(FatalError);


    m1 = Zero;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != Type(Zero))
                    FatalErrorInFunction
                        << "test 2a failed" << abort(FatalError);

    m1 = pTraits<Type>::one*2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != pTraits<Type>::one*2)
                    FatalErrorInFunction
                        << "test 2b failed" << abort(FatalError);

    m1 = m2;
    m1 += m2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 3a failed" << abort(FatalError);

    m1 += (2*m2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 4.0*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 3b failed" << abort(FatalError);

    m1 -= m2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 3.0*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 4a failed" << abort(FatalError);

    m1 -= (2*m2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 4b failed" << abort(FatalError);


    m1 *= s1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 5a failed" << abort(FatalError);

    m1 *= (2*s1);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k)*Foam::sqr(s1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 5b failed" << abort(FatalError);

    m1 /= s1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 6a failed" << abort(FatalError);

    m1 /= (2*s1);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 6b failed" << abort(FatalError);

    m1 += pTraits<Type>::one;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k)+pTraits<Type>::one)
                    FatalErrorInFunction
                        << "test 7a failed" << abort(FatalError);

    m1 -= pTraits<Type>::one;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 7b failed" << abort(FatalError);

    m1 *= scalar(2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 8a failed" << abort(FatalError);

    m1 /= scalar(2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllBlock(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 8b failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testPrimitiveFunctions(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    List<Type> sm(m1[0].size(), pTraits<Type>::zero);

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);

        if (l == 0)
        {
            sm[d] += m1[l][d](i,j,k);
        }
    }

    List<Type> av(sm);

    forAll(m1[0], d)
        av[d] /= cmptProduct(m1[0][d].N());

    m2 = m1;

    s1 = mag(m1)+scalar(1);

    forAll(s1, l)
        forAll(s1[l], d)
            forAllBlock(s1[l][d], i, j, k)
                if
                (
                    s1[l][d](i,j,k)
                 != (
                        Foam::mag(pTraits<Type>::one*(l+d+i+j+k))
                      + scalar(1)
                    )
                )
                    FatalErrorInFunction
                        << "test 9a failed" << abort(FatalError);

    s1 = mag(m1*2)+scalar(1);

    forAll(s1, l)
        forAll(s1[l], d)
            forAllBlock(s1[l][d], i, j, k)
                if
                (
                    s1[l][d](i,j,k)
                 != (
                        Foam::mag(2*pTraits<Type>::one*(l+d+i+j+k))
                      + scalar(1)
                    )
                )
                    FatalErrorInFunction
                        << "test 9b failed" << abort(FatalError);


    List<Type> mx(s1[0].size());

    forAll(m1[0], d)
        mx[d] =
            pTraits<Type>::one*(d + cmptSum(m1[0][d].N()) - 3);

    if (max(m1) != mx)
        FatalErrorInFunction
            << "Test 10a failed" << abort(FatalError);

    if (max(2*m1) != 2*mx)
        FatalErrorInFunction
            << "Test 10b failed" << abort(FatalError);

    if (gMax(m1) != mx)
        FatalErrorInFunction
            << "Test 10c failed" << abort(FatalError);

    if (gMax(2*m1) != 2*mx)
        FatalErrorInFunction
            << "Test 10d failed" << abort(FatalError);


    List<Type> mn(s1[0].size());

    forAll(m1[0], d)
        mn[d] =
            pTraits<Type>::one*d;

    if (min(m1) != mn)
        FatalErrorInFunction
            << "Test 11a failed" << abort(FatalError);

    if (min(-2*m1) != -2*mx)
        FatalErrorInFunction
            << "Test 11b failed" << abort(FatalError);

    if (gMin(m1) != mn)
        FatalErrorInFunction
            << "Test 11c failed" << abort(FatalError);

    if (gMin(-2*m1) != -2*mx)
        FatalErrorInFunction
            << "Test 11d failed" << abort(FatalError);


    if (sum(m1) != sm)
        FatalErrorInFunction
            << "Test 12a failed" << abort(FatalError);

    if (sum(2*m1) != 2*sm)
        FatalErrorInFunction
            << "Test 12b failed" << abort(FatalError);

    if (gSum(m1) != sm*Pstream::nProcs())
        FatalErrorInFunction
            << "Test 12c failed" << abort(FatalError);

    if (gSum(2*m1) != 2*sm*Pstream::nProcs())
        FatalErrorInFunction
            << "Test 12d failed" << abort(FatalError);

    if (average(m1) != av)
        FatalErrorInFunction
            << "Test 13a failed" << abort(FatalError);

    // Multiplying by 2 can give rounding errors for integers

    if (average(1*m1) != av)
        FatalErrorInFunction
            << "Test 13b failed" << abort(FatalError);

    if (gAverage(m1) != av)
        FatalErrorInFunction
            << "Test 13c failed" << abort(FatalError);

    if (gAverage(1*m1) != av)
        FatalErrorInFunction
            << "Test 13d failed" << abort(FatalError);


    sumProd(m1, m1);
    sumProd(m1*2, m1);
    sumProd(m1, m1*2);
    sumProd(m1*2, m1*2);

    meshField<Type,MeshType> m3(max(m1,m2));

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14a failed" << abort(FatalError);

    m3 = max(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14b failed" << abort(FatalError);

    m3 = max(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14c failed" << abort(FatalError);

    m3 = max(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14d failed" << abort(FatalError);


    m3 = min(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15a failed" << abort(FatalError);

    m3 = min(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15b failed" << abort(FatalError);

    m3 = min(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15c failed" << abort(FatalError);

    m3 = min(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15d failed" << abort(FatalError);


    m3 = max(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16a failed" << abort(FatalError);

    m3 = max(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m1[l][d](i,j,k)*2,pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16b failed" << abort(FatalError);

    m3 = max(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m2[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16c failed" << abort(FatalError);

    m3 = max(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m2[l][d](i,j,k)*2,pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16d failed" << abort(FatalError);


    m3 = min(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);


    m3 = -m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != -m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 18a failed" << abort(FatalError);


    m3 = m1*s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19a failed" << abort(FatalError);

    m3 = (m1*2)*s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19b failed" << abort(FatalError);

    m3 = s1*m1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19c failed" << abort(FatalError);

    m3 = s1*(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19d failed" << abort(FatalError);


    m3 = m1/s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != Type(m1[l][d](i,j,k)/s1[l][d](i,j,k)))
                {
                    Pout<< s1[l][d](i,j,k) << endl;
                    FatalErrorInFunction
                        << "test 20a failed" << abort(FatalError);
                }

    m3 = (m1*2)/s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != Type((2*m1[l][d](i,j,k))/s1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 20b failed" << abort(FatalError);


    m3 = m1+m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)+m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21a failed" << abort(FatalError);
    m3 = m1-m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)-m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21b failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testVectorSpaceFunctions(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> m1("m1", fvMsh);
    meshField<Type,MeshType> m2("m2", fvMsh);
    meshField<Type,MeshType> m3("m3", fvMsh);
    meshField<scalar,MeshType> s1("s1", fvMsh);

    label c = 0;

    scalarList m1sm(m1[0].size(), Zero);
    List<Type> m1m1scp(m1[0].size(), Zero);
    List<Type> m1scm(m1[0].size(), Zero);

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;

        if (l == 0)
        {
            m1sm[d] += Foam::mag(m1[l][d](i,j,k));
            m1m1scp[d] += Foam::cmptMultiply(m1[l][d](i,j,k),m1[l][d](i,j,k));
            m1scm[d] += Foam::cmptMag(m1[l][d](i,j,k));
        }
    }

    s1 = cmptMax(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMax(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 22a failed" << abort(FatalError);

    s1 = cmptMax(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMax(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 22a failed" << abort(FatalError);


    s1 = cmptMin(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMin(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 23a failed" << abort(FatalError);

    s1 = cmptMin(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMin(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 23b failed" << abort(FatalError);


    s1 = cmptAv(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptAv(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 24a failed" << abort(FatalError);

    s1 = cmptAv(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptAv(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 24b failed" << abort(FatalError);


    m3 = cmptMag(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMag(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 25a failed" << abort(FatalError);

    m3 = cmptMag(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMag(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 25b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (maxMagSqr(m1)[d] != Foam::magSqr(m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 26a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (maxMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 26b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (gMaxMagSqr(m1)[d] != Foam::magSqr(m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 27a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (gMaxMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 27b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (minMagSqr(m1)[d] != Foam::magSqr(m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 28a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (minMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 28b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (gMinMagSqr(m1)[d] != Foam::magSqr(m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 29a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (gMinMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 29b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (sumMag(m1)[d] != m1sm[d])
            FatalErrorInFunction << "test 30a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (sumMag(2*m1)[d] != 2*m1sm[d])
            FatalErrorInFunction << "test 30b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (gSumMag(m1)[d] != Pstream::nProcs()*m1sm[d])
            FatalErrorInFunction << "test 31a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (gSumMag(2*m1)[d] != Pstream::nProcs()*2*m1sm[d])
            FatalErrorInFunction << "test 31b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (sumCmptProd(m1, m1)[d] != m1m1scp[d])
        {
            Pout<< sumCmptProd(m1, m1)[d] << " " << m1m1scp[d] << endl;
            FatalErrorInFunction << "test 32a failed" << abort(FatalError);
        }

    forAll(m1[0], d)
        if (sumCmptProd(2*m1, m1)[d] != 2*m1m1scp[d])
            FatalErrorInFunction << "test 32b failed" << abort(FatalError);

    forAll(m1[0], d)
        if (sumCmptProd(m1, 2*m1)[d] != 2*m1m1scp[d])
            FatalErrorInFunction << "test 32c failed" << abort(FatalError);

    forAll(m1[0], d)
        if (sumCmptProd(2*m1, 2*m1)[d] != 4*m1m1scp[d])
            FatalErrorInFunction << "test 32d failed" << abort(FatalError);


    forAll(m1[0], d)
        if (sumCmptMag(m1)[d] != m1scm[d])
            FatalErrorInFunction << "test 33a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (sumCmptMag(2*m1)[d] != 2*m1scm[d])
            FatalErrorInFunction << "test 33b failed" << abort(FatalError);

    forAll(m1[0], d)
        if (gSumCmptMag(m1)[d] != Pstream::nProcs()*m1scm[d])
            FatalErrorInFunction << "test 33c failed" << abort(FatalError);

    forAll(m1[0], d)
        if (gSumCmptMag(2*m1)[d] != Pstream::nProcs()*2*m1scm[d])
            FatalErrorInFunction << "test 33d failed" << abort(FatalError);


    m3 = cmptMultiply(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34a failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34b failed" << abort(FatalError);

    m3 = cmptMultiply(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34c failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34d failed" << abort(FatalError);


    m3 = cmptDivide(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35a failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35b failed" << abort(FatalError);

    m3 = cmptDivide(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35c failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35d failed" << abort(FatalError);


    m3 = cmptMultiply(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 36a failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 36b failed" << abort(FatalError);

    m3 = cmptMultiply(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(pTraits<Type>::one,m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 36c failed" << abort(FatalError);

    m3 = cmptMultiply(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(pTraits<Type>::one,2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 36d failed" << abort(FatalError);

    m3 = cmptDivide(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 37a failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 37b failed" << abort(FatalError);

    m3 = cmptDivide(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(pTraits<Type>::one,m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 37c failed" << abort(FatalError);

    m3 = cmptDivide(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllBlock(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(pTraits<Type>::one,2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 37d failed" << abort(FatalError);
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
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1[l][d](i,j,k) = c++;
    }

    -m1;

    m1+m2;
    (2*m1)+m2;
    m1+(2*m2);
    (2*m1)+(2*m2);

    m1-m2;
    (2*m1)-m2;
    m1-(2*m2);
    (2*m1)-(2*m2);

    m1-m2;
    (2*m1)-m2;
    m1-(2*m2);
    (2*m1)-(2*m2);

    m1*s1;
    s1*m1;
    (2*m1)*s1;
    s1*(2*m1);

    m1*(2*s1);
    (2*s1)*m1;
    (2*m1)*(2*s1);
    (2*s1)*(2*m1);

    m1/s1;
    (2*m1)/s1;

    m1/(2*s1);
    (2*m1)/(2*s1);
}

template<class MeshType>
void testScalarFunctions(const fvMesh& fvMsh)
{
    meshField<scalar,MeshType> m1("m1", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = scalar(c+++1);
    }

    m1/m1;
    (m1*2)/m1;
    m1/(m1*2);
}

template<class MeshType>
void testVectorFunctions(const fvMesh& fvMsh)
{
    meshField<vector,MeshType> m1("m1", fvMsh);
    meshField<vector,MeshType> m2("m2", fvMsh);

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<vector>::one*c++;
        m2[l][d](i,j,k) = pTraits<vector>::one*c++;
    }

    m1*m2;
    (m1*2)*m2;
    m1*(m2*2);
    (m1*2)*(m2*2);

    m1 & m2;
    (m1*2) & m2;
    m1 & (m2*2);
    (m1*2) & (m2*2);

    m1 ^ m2;
    (m1*2) ^ m2;
    m1 ^ (m2*2);
    (m1*2) ^ (m2*2);
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
    forAllBlock(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<tensor>::one*c++;
        m2[l][d](i,j,k) = pTraits<symmTensor>::one*c++;
        b3[l][d](i,j,k) = pTraits<sphericalTensor>::one*c++;
        b4[l][d](i,j,k) = pTraits<diagTensor>::one*c++;
        v1[l][d](i,j,k) = vector(l,l+1,l+2); c++;
    }

    m1 & m1;
    (m1*2) & m1;
    m1 & (m1*2);
    (m1*2) & (m1*2);

    m1 & m2;
    (m1*2) & m2;
    m1 & (m2*2);
    (m1*2) & (m2*2);

    m1 && m1;
    (m1*2) && m1;
    m1 && (m1*2);
    (m1*2) && (m1*2);

    m1 && m2;
    (m1*2) && m2;
    m1 && (m2*2);
    (m1*2) && (m2*2);

    m1 & v1;
    (m1*2) & v1;
    v1 & (m1*2);
    (v1*2) & (m1*2);

    m2 & v1;
    (m2*2) & v1;
    v1 & (m2*2);
    (v1*2) & (m2*2);

    b3 & v1;
    (b3*2) & v1;
    v1 & (b3*2);
    (v1*2) & (b3*2);

    b4 & v1;
    (b4*2) & v1;
    v1 & (b4*2);
    (v1*2) & (b4*2);
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
