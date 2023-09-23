#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testConstructors(const fvMesh& fvMsh, const bool deep)
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
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
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
void testIndexing(const fvMesh& fvMsh, const bool deep)
{
    meshField<Type,MeshType> m1
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    label c = 0;

    forAllLevels(m1, l, d, i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
    }

    forAllLevels(m1, l, d, i, j, k)
    {
        if (m1[l][d](i,j,k) != m1[l][d](labelVector(i,j,k)))
            FatalErrorInFunction << "test 1a failed" << abort(FatalError);
    }

    // Direct access on field

    forAllLevels(m1, l, d, i, j, k)
    {
        if (m1[l][d](i,j,k) != m1(l,d,i,j,k))
            FatalErrorInFunction << "test 1b failed" << abort(FatalError);
    }

    forAllLevels(m1, l, d, i, j, k)
    {
        if (m1[l][d](i,j,k) != m1(l,d,labelVector(i,j,k)))
            FatalErrorInFunction << "test 1c failed" << abort(FatalError);
    }

    // Direct access to first level

    forAllDirections(m1, d, i, j, k)
    {
        if (m1[0][d](i,j,k) != m1(d,i,j,k))
            FatalErrorInFunction << "test 1d failed" << abort(FatalError);
    }

    forAllDirections(m1, d, i, j, k)
    {
        if (m1[0][d](i,j,k) != m1(d,labelVector(i,j,k)))
            FatalErrorInFunction << "test 1e failed" << abort(FatalError);
    }

    // Direct access to first direction of first level

    forAllCells(m1, i, j, k)
    {
        if (m1[0][0](i,j,k) != m1(i,j,k))
            FatalErrorInFunction << "test 1f failed" << abort(FatalError);
    }

    forAllCells(m1, i, j, k)
    {
        if (m1[0][0](i,j,k) != m1(labelVector(i,j,k)))
            FatalErrorInFunction << "test 1g failed" << abort(FatalError);
    }
}

template<class Type, class MeshType>
void testMemberOperators(const fvMesh& fvMsh, const bool deep)
{
    meshField<Type,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m1o
    (
        "m1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<Type,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m2o
    (
        "m2o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<scalar,MeshType> s1
    (
        "s1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<scalar,MeshType> s1o
    (
        "s1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);
        s1[l][d](i,j,k) = scalar(l+d+i+j+k+1);
    }

    forAll(m1o, l)
    forAll(m1o[l], d)
    forAllCells(m1o[l][d], i, j, k)
    {
        m1o[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);
        s1o[l][d](i,j,k) = scalar(l+d+i+j+k+1);
    }


    m2 = m1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1a failed" << abort(FatalError);

    m2 = 1.0*m1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1b failed" << abort(FatalError);

    m2 = m1o;
    m2o = m1o;

    forAll(m1o, l)
        forAll(m1o[l], d)
            forAllCells(m1o[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1o[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1c failed" << abort(FatalError);

    m2 = 1.0*m1o;
    m2o = 1.0*m1o;

    forAll(m1o, l)
        forAll(m1o[l], d)
            forAllCells(m1o[l][d], i, j, k)
                if (m2[l][d](i,j,k) != m1o[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 1d failed" << abort(FatalError);

    // Restore
    m2 = m1;


    m1 = Zero;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != Type(Zero))
                    FatalErrorInFunction
                        << "test 2a failed" << abort(FatalError);

    m1 = pTraits<Type>::one*2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != pTraits<Type>::one*2)
                    FatalErrorInFunction
                        << "test 2b failed" << abort(FatalError);

    m1 = m2;
    m1 += m2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 3a failed" << abort(FatalError);

    m1 += (2*m2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 4.0*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 3b failed" << abort(FatalError);

    if (deep)
    {
        m1o = m2o;
        m1o += m2;

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 3c failed" << abort(FatalError);

        m1o += (2*m2);

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != 4.0*m2[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 3d failed" << abort(FatalError);
    }

    m1 -= m2;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 3.0*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 4a failed" << abort(FatalError);

    m1 -= (2*m2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 4b failed" << abort(FatalError);

    if (deep)
    {
        m1o -= m2;

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != 3.0*m2[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 4c failed" << abort(FatalError);

        m1o -= (2*m2);

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != m2[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 4d failed" << abort(FatalError);
    }

    m1 *= s1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 5a failed" << abort(FatalError);

    m1 *= (2*s1);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k)*Foam::sqr(s1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 5b failed" << abort(FatalError);

    if (deep)
    {
        m1o *= s1;

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != m2[l][d](i,j,k)*s1[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 5c failed" << abort(FatalError);

        m1o *= (2*s1);

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != 2*m2[l][d](i,j,k)*Foam::sqr(s1[l][d](i,j,k)))
                        FatalErrorInFunction
                            << "test 5d failed" << abort(FatalError);
    }

    m1 /= s1;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 6a failed" << abort(FatalError);

    m1 /= (2*s1);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 6b failed" << abort(FatalError);

    if (deep)
    {
        m1o /= s1;

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != 2*m2[l][d](i,j,k)*s1[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 6c failed" << abort(FatalError);

        m1o /= (2*s1);

        forAll(m1o, l)
            forAll(m1o[l], d)
                forAllCells(m1o[l][d], i, j, k)
                    if (m1o[l][d](i,j,k) != m2[l][d](i,j,k))
                        FatalErrorInFunction
                            << "test 6d failed" << abort(FatalError);
    }

    m1 += pTraits<Type>::one;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k)+pTraits<Type>::one)
                    FatalErrorInFunction
                        << "test 7a failed" << abort(FatalError);

    m1 -= pTraits<Type>::one;

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 7b failed" << abort(FatalError);

    m1 *= scalar(2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 8a failed" << abort(FatalError);

    m1 /= scalar(2);

    forAll(m1, l)
        forAll(m1[l], d)
            forAllCells(m1[l][d], i, j, k)
                if (m1[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 8b failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testPrimitiveFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<Type,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m1o
    (
        "m1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<Type,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m2o
    (
        "m2o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<scalar,MeshType> s1
    (
        "s1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<scalar,MeshType> s1o
    (
        "s1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    List<Type> sm(m1[0].size(), pTraits<Type>::zero);

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);

        if (l == 0)
        {
            sm[d] += m1[l][d](i,j,k);
        }
    }

    forAll(m1o, l)
    forAll(m1o[l], d)
    forAllCells(m1o[l][d], i, j, k)
    {
        m1o[l][d](i,j,k) = pTraits<Type>::one*(l+d+i+j+k);
    }

    List<Type> av(sm);

    forAll(m1[0], d)
        av[d] /= cmptProduct(m1[0][d].N());

    m2 = m1;
    m2o = m1o;

    s1 = mag(m1)+scalar(1);
    s1o = mag(m1o)+scalar(1);

    forAll(s1, l)
        forAll(s1[l], d)
            forAllCells(s1[l][d], i, j, k)
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
    s1o = mag(m1o*2)+scalar(1);

    forAll(s1, l)
        forAll(s1[l], d)
            forAllCells(s1[l][d], i, j, k)
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
    {
        mx[d] =
            pTraits<Type>::one*(d + cmptSum(m1[0][d].N()) - 3);
    }

    List<Type> gmx(mx);

    forAll(m1[0], d)
        reduce(gmx[d], maxOp<Type>());

    if (max(m1) != mx)
        FatalErrorInFunction
            << "Test 10a failed" << abort(FatalError);

    if (max(2*m1) != 2*mx)
        FatalErrorInFunction
            << "Test 10b failed" << abort(FatalError);

    if (gMax(m1) != gmx)
        FatalErrorInFunction
            << "Test 10c failed" << abort(FatalError);

    if (gMax(2*m1) != 2*gmx)
        FatalErrorInFunction
            << "Test 10d failed" << abort(FatalError);


    List<Type> mn(s1[0].size());

    forAll(m1[0], d)
    {
        mn[d] =
            pTraits<Type>::one*d;
    }

    List<Type> gmn(mn);

    forAll(m1[0], d)
        reduce(gmn[d], maxOp<Type>());

    if (min(m1) != mn)
        FatalErrorInFunction
            << "Test 11a failed" << abort(FatalError);

    if (min(-2*m1) != -2*mx)
        FatalErrorInFunction
            << "Test 11b failed" << abort(FatalError);

    if (gMin(m1) != gmn)
        FatalErrorInFunction
            << "Test 11c failed" << abort(FatalError);

    if (gMin(-2*m1) != -2*gmx)
        FatalErrorInFunction
            << "Test 11d failed" << abort(FatalError);


    List<Type> gsm(sm);

    forAll(gsm, d)
        reduce(gsm[d], sumOp<Type>());

    if (sum(m1) != sm)
        FatalErrorInFunction
            << "Test 12a failed" << abort(FatalError);

    if (sum(2*m1) != 2*sm)
        FatalErrorInFunction
            << "Test 12b failed" << abort(FatalError);

    if (gSum(m1) != gsm)
        FatalErrorInFunction
            << "Test 12c failed" << abort(FatalError);

    if (gSum(2*m1) != 2*gsm)
        FatalErrorInFunction
            << "Test 12d failed" << abort(FatalError);


    List<Type> gav(gsm);

    forAll(gav, d)
    {
        scalar N = cmptProduct(m1[0][d].N());
        reduce(N, sumOp<scalar>());

        gav[d] /= N;
    }

    if (average(m1) != av)
        FatalErrorInFunction
            << "Test 13a failed" << abort(FatalError);

    // Note: multiplying by 2 can give rounding errors for integers

    if (average(1*m1) != av)
        FatalErrorInFunction
            << "Test 13b failed" << abort(FatalError);

    if (gAverage(m1) != gav)
        FatalErrorInFunction
            << "Test 13c failed" << endl
            << gAverage(m1) << endl
            << gav << endl << abort(FatalError);

    if (gAverage(1*m1) != gav)
        FatalErrorInFunction
            << "Test 13d failed" << abort(FatalError);


    sumProd(m1, m1);
    sumProd(m1*2, m1);
    sumProd(m1, m1*2);
    sumProd(m1*2, m1*2);

    sumProd(m1, m1o);
    sumProd(m1*2, m1o);
    sumProd(m1, m1o*2);
    sumProd(m1*2, m1o*2);

    meshField<Type,MeshType> m3(max(m1,m2));

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14a failed" << abort(FatalError);

    m3 = max(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14b failed" << abort(FatalError);

    m3 = max(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14c failed" << abort(FatalError);

    m3 = max(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14d failed" << abort(FatalError);

    meshField<Type,MeshType> m3o(max(m1,m2o));

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14e failed" << abort(FatalError);

    m3o = max(m1*2,m2o);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14f failed" << abort(FatalError);

    m3o = max(m1,m2o*2);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != 2*m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14g failed" << abort(FatalError);

    m3o = max(m1*2,m2o*2);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 14h failed" << abort(FatalError);


    m3 = min(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15a failed" << abort(FatalError);

    m3 = min(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15b failed" << abort(FatalError);

    m3 = min(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15c failed" << abort(FatalError);

    m3 = min(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15d failed" << abort(FatalError);

    m3o = min(m1,m2o);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15e failed" << abort(FatalError);

    m3o = min(m1*2,m2o);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15f failed" << abort(FatalError);

    m3o = min(m1,m2o*2);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15g failed" << abort(FatalError);

    m3o = min(m1*2,m2o*2);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != 2*m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 15h failed" << abort(FatalError);


    m3 = max(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16a failed" << abort(FatalError);

    m3 = max(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m1[l][d](i,j,k)*2,pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16b failed" << abort(FatalError);

    m3 = max(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m2[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16c failed" << abort(FatalError);

    m3 = max(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != max(m2[l][d](i,j,k)*2,pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 16d failed" << abort(FatalError);


    m3 = min(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);

    m3 = min(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != min(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 17a failed" << abort(FatalError);


    m3 = -m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != -m1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 18a failed" << abort(FatalError);


    m3 = m1*s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19a failed" << abort(FatalError);

    m3 = (m1*2)*s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19b failed" << abort(FatalError);

    m3 = s1*m1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19c failed" << abort(FatalError);

    m3 = s1*(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19d failed" << abort(FatalError);

    m3o = m1*s1o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19e failed" << abort(FatalError);

    m3o = (m1*2)*s1o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19f failed" << abort(FatalError);

    m3o = s1*m1o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 19g failed" << abort(FatalError);

    m3o = s1*(m1o*2);

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)*s1[l][d](i,j,k)*2)
                    FatalErrorInFunction
                        << "test 19h failed" << abort(FatalError);

    m3 = m1/s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != Type(m1[l][d](i,j,k)/s1[l][d](i,j,k)))
                {
                    Pout<< s1[l][d](i,j,k) << endl;
                    FatalErrorInFunction
                        << "test 20a failed" << abort(FatalError);
                }

    m3 = (m1*2)/s1;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != Type((2*m1[l][d](i,j,k))/s1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 20b failed" << abort(FatalError);

    m3o = m1/s1o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != Type(m1[l][d](i,j,k)/s1[l][d](i,j,k)))
                {
                    Pout<< s1[l][d](i,j,k) << endl;
                    FatalErrorInFunction
                        << "test 20c failed" << abort(FatalError);
                }

    m3o = (m1*2)/s1o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != Type((2*m1[l][d](i,j,k))/s1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 20d failed" << abort(FatalError);

    m3 = m1+m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)+m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21a failed" << abort(FatalError);
    m3 = m1-m2;

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != m1[l][d](i,j,k)-m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21b failed" << abort(FatalError);

    m3o = m1+m2o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)+m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21c failed" << abort(FatalError);
    m3o = m1-m2o;

    forAll(m3o, l)
        forAll(m3o[l], d)
            forAllCells(m3o[l][d], i, j, k)
                if (m3o[l][d](i,j,k) != m1[l][d](i,j,k)-m2[l][d](i,j,k))
                    FatalErrorInFunction
                        << "test 21d failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testVectorSpaceFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<Type,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m3
    (
        "m3",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<scalar,MeshType> s1
    (
        "s1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    label c = 0;

    scalarList m1sm(m1[0].size(), Zero);
    List<Type> m1m1scp(m1[0].size(), Zero);
    List<Type> m1scm(m1[0].size(), Zero);

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
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
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMax(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 22a failed" << abort(FatalError);

    s1 = cmptMax(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMax(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 22a failed" << abort(FatalError);


    s1 = cmptMin(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMin(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 23a failed" << abort(FatalError);

    s1 = cmptMin(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptMin(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 23b failed" << abort(FatalError);


    s1 = cmptAv(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptAv(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 24a failed" << abort(FatalError);

    s1 = cmptAv(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (s1[l][d](i,j,k) != cmptAv(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 24b failed" << abort(FatalError);


    m3 = cmptMag(m1);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMag(m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 25a failed" << abort(FatalError);

    m3 = cmptMag(m1*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMag(2*m1[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 25b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (maxMagSqr(m1)[d] != Foam::magSqr(m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 26a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (maxMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 26b failed" << abort(FatalError);


    // For this and subsequent global operators, only perform for colocated
    // meshes which is easier to verify. If it works on colocated meshes, it
    // should also work on other meshes...

    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gMaxMagSqr(m1)[d] != Foam::magSqr(m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 27a failed" << abort(FatalError);

    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gMaxMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](m1[0][d].N()-unitXYZ)))
            FatalErrorInFunction << "test 27b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (minMagSqr(m1)[d] != Foam::magSqr(m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 28a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (minMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 28b failed" << abort(FatalError);


    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gMinMagSqr(m1)[d] != Foam::magSqr(m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 29a failed" << abort(FatalError);

    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gMinMagSqr(2*m1)[d] != Foam::magSqr(2*m1[0][d](0,0,0)))
            FatalErrorInFunction << "test 29b failed" << abort(FatalError);


    forAll(m1[0], d)
        if (sumMag(m1)[d] != m1sm[d])
            FatalErrorInFunction << "test 30a failed" << abort(FatalError);

    forAll(m1[0], d)
        if (sumMag(2*m1)[d] != 2*m1sm[d])
            FatalErrorInFunction << "test 30b failed" << abort(FatalError);


    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gSumMag(m1)[d] != Pstream::nProcs()*m1sm[d])
            FatalErrorInFunction << "test 31a failed" << abort(FatalError);

    if (MeshType::numberOfDirections == 1)
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

    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gSumCmptMag(m1)[d] != Pstream::nProcs()*m1scm[d])
            FatalErrorInFunction << "test 33c failed" << abort(FatalError);

    if (MeshType::numberOfDirections == 1)
    forAll(m1[0], d)
        if (gSumCmptMag(2*m1)[d] != Pstream::nProcs()*2*m1scm[d])
            FatalErrorInFunction << "test 33d failed" << abort(FatalError);


    m3 = cmptMultiply(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34a failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34b failed" << abort(FatalError);

    m3 = cmptMultiply(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34c failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 34d failed" << abort(FatalError);


    m3 = cmptDivide(m1,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35a failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35b failed" << abort(FatalError);

    m3 = cmptDivide(m1,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35c failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 35d failed" << abort(FatalError);


    m3 = cmptMultiply(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 36a failed" << abort(FatalError);

    m3 = cmptMultiply(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 36b failed" << abort(FatalError);

    m3 = cmptMultiply(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(pTraits<Type>::one,m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 36c failed" << abort(FatalError);

    m3 = cmptMultiply(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptMultiply(pTraits<Type>::one,2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 36d failed" << abort(FatalError);

    m3 = cmptDivide(m1,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 37a failed" << abort(FatalError);

    m3 = cmptDivide(m1*2,pTraits<Type>::one);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(2*m1[l][d](i,j,k),pTraits<Type>::one))
                    FatalErrorInFunction
                        << "test 37b failed" << abort(FatalError);

    m3 = cmptDivide(pTraits<Type>::one,m2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(pTraits<Type>::one,m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 37c failed" << abort(FatalError);

    m3 = cmptDivide(pTraits<Type>::one,m2*2);

    forAll(m3, l)
        forAll(m3[l], d)
            forAllCells(m3[l][d], i, j, k)
                if (m3[l][d](i,j,k) != cmptDivide(pTraits<Type>::one,2*m2[l][d](i,j,k)))
                    FatalErrorInFunction
                        << "test 37d failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testStencilFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<Type,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<scalar,MeshType> s1
    (
        "s1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<Type,MeshType> m1o
    (
        "m1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<Type,MeshType> m2o
    (
        "m2o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    meshField<scalar,MeshType> s1o
    (
        "s1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1[l][d](i,j,k) = c++;
    }

    forAll(m1o, l)
    forAll(m1o[l], d)
    forAllCells(m1o[l][d], i, j, k)
    {
        m1o[l][d](i,j,k) = pTraits<Type>::one*c++;
        m2o[l][d](i,j,k) = pTraits<Type>::one*c++;
        s1o[l][d](i,j,k) = c++;
    }

    -m1;

    m1+m2;
    (2*m1)+m2;
    m1+(2*m2);
    (2*m1)+(2*m2);

    m1+m2o;
    (2*m1)+m2o;
    m1+(2*m2o);
    (2*m1)+(2*m2o);

    m1-m2;
    (2*m1)-m2;
    m1-(2*m2);
    (2*m1)-(2*m2);

    m1-m2o;
    (2*m1)-m2o;
    m1-(2*m2o);
    (2*m1)-(2*m2o);

    m1-m2;
    (2*m1)-m2;
    m1-(2*m2);
    (2*m1)-(2*m2);

    m1-m2o;
    (2*m1)-m2o;
    m1-(2*m2o);
    (2*m1)-(2*m2o);

    m1*s1;
    s1*m1;
    (2*m1)*s1;
    s1*(2*m1);

    m1*s1o;
    s1*m1o;
    (2*m1)*s1o;
    s1*(2*m1o);

    m1*(2*s1);
    (2*s1)*m1;
    (2*m1)*(2*s1);
    (2*s1)*(2*m1);

    m1*(2*s1o);
    (2*s1)*m1o;
    (2*m1)*(2*s1o);
    (2*s1)*(2*m1o);

    m1/s1;
    (2*m1)/s1;

    m1/s1o;
    (2*m1)/s1o;

    m1/(2*s1);
    (2*m1)/(2*s1);

    m1/(2*s1o);
    (2*m1)/(2*s1o);
}

template<class MeshType>
void testScalarFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<scalar,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<scalar,MeshType> m1o
    (
        "m1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        !deep
    );

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
    {
        m1[l][d](i,j,k) = scalar(c+++1);
    }

    forAll(m1o, l)
    forAll(m1o[l], d)
    forAllCells(m1o[l][d], i, j, k)
    {
        m1o[l][d](i,j,k) = scalar(c+++1);
    }

    m1/m1;
    (m1*2)/m1;
    m1/(m1*2);

    m1/m1o;
    (m1*2)/m1o;
    m1/(m1o*2);
}

template<class MeshType>
void testVectorFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<vector,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<vector,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    label c = 0;

    forAll(m1, l)
    forAll(m1[l], d)
    forAllCells(m1[l][d], i, j, k)
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
void testTensorFunctions(const fvMesh& fvMsh, const bool deep)
{
    meshField<tensor,MeshType> m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<symmTensor,MeshType> m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<sphericalTensor,MeshType> b3
    (
        "b3",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<diagTensor,MeshType> b4
    (
        "b4",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

    meshField<vector,MeshType> v1
    (
        "v1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        deep
    );

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

    for (int deep = 0; deep < 2; deep++)
    {
        testConstructors<label,colocated>(fvMsh, deep);
        testConstructors<scalar,colocated>(fvMsh, deep);

        testConstructors<vector,colocated>(fvMsh, deep);
        testConstructors<tensor,colocated>(fvMsh, deep);
        testConstructors<symmTensor,colocated>(fvMsh, deep);
        testConstructors<sphericalTensor,colocated>(fvMsh, deep);
        testConstructors<diagTensor,colocated>(fvMsh, deep);

        testConstructors<faceScalar,colocated>(fvMsh, deep);
        testConstructors<edgeScalar,colocated>(fvMsh, deep);
        testConstructors<vertexScalar,colocated>(fvMsh, deep);

        testConstructors<faceVector,colocated>(fvMsh, deep);
        testConstructors<edgeVector,colocated>(fvMsh, deep);
        testConstructors<vertexVector,colocated>(fvMsh, deep);

        testConstructors<stencil,colocated>(fvMsh, deep);
        testConstructors<diagStencil,colocated>(fvMsh, deep);


        testIndexing<label,colocated>(fvMsh, deep);
        testIndexing<scalar,colocated>(fvMsh, deep);

        testIndexing<vector,colocated>(fvMsh, deep);
        testIndexing<tensor,colocated>(fvMsh, deep);
        testIndexing<symmTensor,colocated>(fvMsh, deep);
        testIndexing<sphericalTensor,colocated>(fvMsh, deep);
        testIndexing<diagTensor,colocated>(fvMsh, deep);

        testIndexing<faceScalar,colocated>(fvMsh, deep);
        testIndexing<edgeScalar,colocated>(fvMsh, deep);
        testIndexing<vertexScalar,colocated>(fvMsh, deep);

        testIndexing<faceVector,colocated>(fvMsh, deep);
        testIndexing<edgeVector,colocated>(fvMsh, deep);
        testIndexing<vertexVector,colocated>(fvMsh, deep);

        testIndexing<stencil,colocated>(fvMsh, deep);
        testIndexing<diagStencil,colocated>(fvMsh, deep);


        testMemberOperators<label,colocated>(fvMsh, deep);
        testMemberOperators<scalar,colocated>(fvMsh, deep);

        testMemberOperators<vector,colocated>(fvMsh, deep);
        testMemberOperators<tensor,colocated>(fvMsh, deep);
        testMemberOperators<symmTensor,colocated>(fvMsh, deep);
        testMemberOperators<sphericalTensor,colocated>(fvMsh, deep);
        testMemberOperators<diagTensor,colocated>(fvMsh, deep);

        testMemberOperators<faceScalar,colocated>(fvMsh, deep);
        testMemberOperators<edgeScalar,colocated>(fvMsh, deep);
        testMemberOperators<vertexScalar,colocated>(fvMsh, deep);

        testMemberOperators<faceVector,colocated>(fvMsh, deep);
        testMemberOperators<edgeVector,colocated>(fvMsh, deep);
        testMemberOperators<vertexVector,colocated>(fvMsh, deep);

        testMemberOperators<stencil,colocated>(fvMsh, deep);
        testMemberOperators<diagStencil,colocated>(fvMsh, deep);


        testPrimitiveFunctions<label,colocated>(fvMsh, deep);
        testPrimitiveFunctions<scalar,colocated>(fvMsh, deep);

        testPrimitiveFunctions<vector,colocated>(fvMsh, deep);
        testPrimitiveFunctions<tensor,colocated>(fvMsh, deep);
        testPrimitiveFunctions<symmTensor,colocated>(fvMsh, deep);
        testPrimitiveFunctions<sphericalTensor,colocated>(fvMsh, deep);
        testPrimitiveFunctions<diagTensor,colocated>(fvMsh, deep);

        testVectorSpaceFunctions<vector,colocated>(fvMsh, deep);
        testVectorSpaceFunctions<tensor,colocated>(fvMsh, deep);
        testVectorSpaceFunctions<symmTensor,colocated>(fvMsh, deep);
        testVectorSpaceFunctions<sphericalTensor,colocated>(fvMsh, deep);
        testVectorSpaceFunctions<diagTensor,colocated>(fvMsh, deep);

        testStencilFunctions<stencil,colocated>(fvMsh, deep);
        testStencilFunctions<diagStencil,colocated>(fvMsh, deep);

        testScalarFunctions<colocated>(fvMsh, deep);
        testVectorFunctions<colocated>(fvMsh, deep);
        testTensorFunctions<colocated>(fvMsh, deep);
    }

    // Staggered

    if (fvMsh.structured())
    {
        for (int deep = 0; deep < 2; deep++)
        {
            testConstructors<label,staggered>(fvMsh, deep);
            testConstructors<scalar,staggered>(fvMsh, deep);
            testConstructors<vector,staggered>(fvMsh, deep);

            testConstructors<tensor,staggered>(fvMsh, deep);
            testConstructors<symmTensor,staggered>(fvMsh, deep);
            testConstructors<sphericalTensor,staggered>(fvMsh, deep);
            testConstructors<diagTensor,staggered>(fvMsh, deep);

            testConstructors<faceScalar,staggered>(fvMsh, deep);
            testConstructors<edgeScalar,staggered>(fvMsh, deep);
            testConstructors<vertexScalar,staggered>(fvMsh, deep);

            testConstructors<faceVector,staggered>(fvMsh, deep);
            testConstructors<edgeVector,staggered>(fvMsh, deep);
            testConstructors<vertexVector,staggered>(fvMsh, deep);

            testConstructors<stencil,staggered>(fvMsh, deep);
            testConstructors<diagStencil,staggered>(fvMsh, deep);


            testIndexing<label,staggered>(fvMsh, deep);
            testIndexing<scalar,staggered>(fvMsh, deep);

            testIndexing<vector,staggered>(fvMsh, deep);
            testIndexing<tensor,staggered>(fvMsh, deep);
            testIndexing<symmTensor,staggered>(fvMsh, deep);
            testIndexing<sphericalTensor,staggered>(fvMsh, deep);
            testIndexing<diagTensor,staggered>(fvMsh, deep);

            testIndexing<faceScalar,staggered>(fvMsh, deep);
            testIndexing<edgeScalar,staggered>(fvMsh, deep);
            testIndexing<vertexScalar,staggered>(fvMsh, deep);

            testIndexing<faceVector,staggered>(fvMsh, deep);
            testIndexing<edgeVector,staggered>(fvMsh, deep);
            testIndexing<vertexVector,staggered>(fvMsh, deep);

            testIndexing<stencil,staggered>(fvMsh, deep);
            testIndexing<diagStencil,staggered>(fvMsh, deep);


            testMemberOperators<label,staggered>(fvMsh, deep);
            testMemberOperators<scalar,staggered>(fvMsh, deep);

            testMemberOperators<vector,staggered>(fvMsh, deep);
            testMemberOperators<tensor,staggered>(fvMsh, deep);
            testMemberOperators<symmTensor,staggered>(fvMsh, deep);
            testMemberOperators<sphericalTensor,staggered>(fvMsh, deep);
            testMemberOperators<diagTensor,staggered>(fvMsh, deep);

            testMemberOperators<faceScalar,staggered>(fvMsh, deep);
            testMemberOperators<edgeScalar,staggered>(fvMsh, deep);
            testMemberOperators<vertexScalar,staggered>(fvMsh, deep);

            testMemberOperators<faceVector,staggered>(fvMsh, deep);
            testMemberOperators<edgeVector,staggered>(fvMsh, deep);
            testMemberOperators<vertexVector,staggered>(fvMsh, deep);

            testMemberOperators<stencil,staggered>(fvMsh, deep);
            testMemberOperators<diagStencil,staggered>(fvMsh, deep);


            testPrimitiveFunctions<label,staggered>(fvMsh, deep);
            testPrimitiveFunctions<scalar,staggered>(fvMsh, deep);

            testPrimitiveFunctions<vector,staggered>(fvMsh, deep);
            testPrimitiveFunctions<tensor,staggered>(fvMsh, deep);
            testPrimitiveFunctions<symmTensor,staggered>(fvMsh, deep);
            testPrimitiveFunctions<sphericalTensor,staggered>(fvMsh, deep);
            testPrimitiveFunctions<diagTensor,staggered>(fvMsh, deep);


            testVectorSpaceFunctions<vector,staggered>(fvMsh, deep);
            testVectorSpaceFunctions<tensor,staggered>(fvMsh, deep);
            testVectorSpaceFunctions<symmTensor,staggered>(fvMsh, deep);
            testVectorSpaceFunctions<sphericalTensor,staggered>(fvMsh, deep);
            testVectorSpaceFunctions<diagTensor,staggered>(fvMsh, deep);


            testStencilFunctions<stencil,staggered>(fvMsh, deep);
            testStencilFunctions<diagStencil,staggered>(fvMsh, deep);


            testScalarFunctions<staggered>(fvMsh, deep);
            testVectorFunctions<staggered>(fvMsh, deep);
            testTensorFunctions<staggered>(fvMsh, deep);
        }
    }
}
