#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "faceFields.H"
#include "stencilFields.H"
#include "diagStencilFields.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SetType>
void testConstructors(const fvMesh& fvMsh, const bool deep)
{
    typedef typename SetType::dataType Type;
    typedef typename SetType::meshType MeshType;

    // Field from name and mesh

    SetType m1("m1", fvMsh);

    // Read

    SetType m2
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        deep
    );

    // Copy from m2 with same name

    SetType m3a(m2);
    SetType m3b(m2, true);

    // Copy from m2 with new name

    SetType m4a("m4a", m2);
    SetType m4b("m4b", m2, true);

    // Copy from tmp of m2 with same name

    SetType m5a(2.0*m2);
    SetType m5b(2.0*m2, true);

    // Copy from tmp m2 with new name

    SetType m6a("m6a", 2.0*m2);
    SetType m6b("m6b", 2.0*m2, true);
}

template<class SetType>
void testMemberOperators(const fvMesh& fvMsh, const bool deep)
{
    SetType m1
    (
        "m1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        deep
    );

    m1 = Zero;

    SetType m1o
    (
        "m1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        !deep
    );

    m1o = Zero;

    SetType m2
    (
        "m2",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        deep
    );

    m2 = Zero;

    SetType m2o
    (
        "m2o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        !deep
    );

    m2o = Zero;

    typedef typename SetType::dataType Type;
    typedef typename SetType::meshType MeshType;
    const direction N = SetType::nSetComponents;

    meshFieldSet<scalar,MeshType,N> s1
    (
        "s1",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        deep
    );

    s1 = 2.0;

    meshFieldSet<scalar,MeshType,N> s1o
    (
        "s1o",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        !deep
    );

    s1o = 2.0;

    m2*2.0;

    m2 = m1;
    m2 = 1.0*m1;

    m2 = m1o;
    m2o = m1o;

    m2 = 1.0*m1o;
    m2o = 1.0*m1o;

    m1 = Zero;
    m1 = pTraits<Type>::one*2;

    m1 = m2;
    m1 += m2;

    m1 += (2*m2);

    if (deep)
    {
        m1o = m2o;
        m1o += m2;
        m1o += (2*m2);
    }

    m1 -= m2;
    m1 -= (2*m2);

    if (deep)
    {
        m1o -= m2;
        m1o -= (2*m2);
    }

    m1 *= s1;
    m1 *= (2*s1);

    if (deep)
    {
        m1o *= s1;
        m1o *= (2*s1);
    }

    m1 /= s1;
    m1 /= (2*s1);

    if (deep)
    {
        m1o /= s1;
        m1o /= (2*s1);
    }

    m1 += pTraits<Type>::one;
    m1 -= pTraits<Type>::one;

    m1 *= scalar(2);
    m1 /= scalar(2);

    -m1;
    -(2.0*m1);

    m1+m2;
    m1+(2.0*m2);
    (2.0*m1)+m2;
    (2.0*m1)+(2.0*m2);

    m1-m2;
    m1-(2.0*m2);
    (2.0*m1)-m2;
    (2.0*m1)-(2.0*m2);

    m1*s1;
    m1*(2.0*s1);
    (2.0*m1)*s1;
    (2.0*m1)*(2.0*s1);

    m1/s1;
    m1/(2.0*s1);
    (2.0*m1)/s1;
    (2.0*m1)/(2.0*s1);
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
        testConstructors<colocatedScalarFaceField>(fvMsh, deep);
        testConstructors<colocatedVectorFaceField>(fvMsh, deep);

        testMemberOperators<colocatedScalarFaceField>(fvMsh, deep);
        testMemberOperators<colocatedVectorFaceField>(fvMsh, deep);


        testConstructors<colocatedScalarStencilField>(fvMsh, deep);
        testConstructors<colocatedVectorStencilField>(fvMsh, deep);

        testMemberOperators<colocatedScalarStencilField>(fvMsh, deep);
        testMemberOperators<colocatedVectorStencilField>(fvMsh, deep);


        testConstructors<colocatedScalarDiagStencilField>(fvMsh, deep);
        testConstructors<colocatedVectorDiagStencilField>(fvMsh, deep);

        testMemberOperators<colocatedScalarDiagStencilField>(fvMsh, deep);
        testMemberOperators<colocatedVectorDiagStencilField>(fvMsh, deep);
    }

    // Staggered

    if (fvMsh.structured())
    {
        for (int deep = 0; deep < 2; deep++)
        {
            testConstructors<staggeredScalarFaceField>(fvMsh, deep);
            testConstructors<staggeredVectorFaceField>(fvMsh, deep);

            testMemberOperators<staggeredScalarFaceField>(fvMsh, deep);
            testMemberOperators<staggeredVectorFaceField>(fvMsh, deep);


            testConstructors<staggeredScalarStencilField>(fvMsh, deep);
            testConstructors<staggeredVectorStencilField>(fvMsh, deep);

            testMemberOperators<staggeredScalarStencilField>(fvMsh, deep);
            testMemberOperators<staggeredVectorStencilField>(fvMsh, deep);


            testConstructors<staggeredScalarDiagStencilField>(fvMsh, deep);
            testConstructors<staggeredVectorDiagStencilField>(fvMsh, deep);

            testMemberOperators<staggeredScalarDiagStencilField>(fvMsh, deep);
            testMemberOperators<staggeredVectorDiagStencilField>(fvMsh, deep);
        }
    }
}
