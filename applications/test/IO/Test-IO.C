#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "IO.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> createField
(
    const word name,
    const fvMesh& fvMsh,
    const bool deep
)
{
    tmp<meshField<Type,MeshType>> tf
    (
        new meshField<Type,MeshType>
        (
            name,
            fvMsh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true,
            deep
        )
    );

    meshField<Type,MeshType>& f = tf.ref();

    f = Zero;

    forAll(f, l)
    forAll(f[l], d)
    {
        const labelVector S = f.I(l,d).lower()-unitXYZ;
        const labelVector E = f.I(l,d).upper()+unitXYZ;

        for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
        for (int k = S.z(); k < E.z(); k++)
        {
            f(l,d,i,j,k) = (Pstream::myProcNo()+i+j+k+l+d)*pTraits<Type>::one;
        }
    }

    return tf;
}

inline bool structured
(
    const label i,
    const label j,
    const label k,
    const labelVector N
)
{
    return
        (i > 0 && i < N.x()-1 && j > 0 && j < N.y()-1)
     || (i > 0 && i < N.x()-1 && k > 0 && k < N.z()-1)
     || (j > 0 && j < N.y()-1 && k > 0 && k < N.z()-1);
}

template<class Type, class MeshType>
void testField(const meshField<Type,MeshType>& f, bool ghosts)
{
    forAll(f, l)
    forAll(f[l], d)
    {
        const labelVector S = f.I(l,d).lower() - unitXYZ*label(ghosts);
        const labelVector E = f.I(l,d).upper() + unitXYZ*label(ghosts);
        const labelVector N(E-S);

        // In the case of an unstructured mesh with ghosts, don't check edge and
        // vertex cells

        for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
        for (int k = S.z(); k < E.z(); k++)
        if
        (
            !ghosts
          || f.fvMsh().structured()
          || structured(i-S.x(), j-S.y(), k-S.z(), N)
        )
        {
            if
            (
                f(l,d,i,j,k)
             != (Pstream::myProcNo()+i+j+k+l+d)*pTraits<Type>::one
            )
            {
                FatalErrorInFunction
                    << "Test failed for " << nl
                    << "ijk = " << labelVector(i,j,k) << nl
                    << "l = " << l << nl
                    << "d = " << d << nl
                    << "type = " << pTraits<Type>::typeName << nl
                    << "meshType = " << MeshType::typeName << nl
                    << "test ghosts = " << Switch(ghosts) << nl
                    << "expected value = "
                    << (Pstream::myProcNo()+i+j+k)*pTraits<Type>::one << nl
                    << "actual value = " << f(l,d,i,j,k)
                    << abort(FatalError);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaIO.H"

    // Macros

    #define CREATEFIELD(TYPE, TYPENAME, MESHTYPE, DEEP)                 \
        tmp<meshField<TYPE,MESHTYPE>> MESHTYPE##TYPENAME##TestField     \
        (                                                               \
            createField<TYPE,MESHTYPE>                                  \
            (                                                           \
                #MESHTYPE #TYPENAME "TestField",                        \
                fvMsh,                                                  \
                DEEP                                                    \
            )                                                           \
        );

    #define RESETFIELD(TYPENAME, MESHTYPE)                              \
        MESHTYPE##TYPENAME##TestField.ref() = Zero;                     \
        MESHTYPE##TYPENAME##TestField.ref().readOpt() =                 \
            IOobject::MUST_READ;

    #define TESTFIELD(TYPE, TYPENAME, MESHTYPE, GHOSTS)                 \
        testField<TYPE,MESHTYPE>                                        \
        (                                                               \
            MESHTYPE##TYPENAME##TestField.ref(),                        \
            GHOSTS                                                      \
        );


    #define TEST(MESHTYPE,TIME,GHOSTS,PARTITIONED,DEEP)                 \
    {                                                                   \
        runTime.setTime(TIME, TIME);                                    \
                                                                        \
        io.ghosts() = GHOSTS;                                           \
        io.partitioned() = PARTITIONED;                                 \
                                                                        \
        CREATEFIELD(label,Label,MESHTYPE,DEEP)                          \
        CREATEFIELD(scalar,Scalar,MESHTYPE,DEEP)                        \
        CREATEFIELD(vector,Vector,MESHTYPE,DEEP)                        \
        CREATEFIELD(tensor,Tensor,MESHTYPE,DEEP)                        \
        CREATEFIELD(diagTensor,DiagTensor,MESHTYPE,DEEP)                \
        CREATEFIELD(symmTensor,SymmTensor,MESHTYPE,DEEP)                \
        CREATEFIELD(sphericalTensor,SphericalTensor,MESHTYPE,DEEP)      \
        CREATEFIELD(faceScalar,FaceScalar,MESHTYPE,DEEP)                \
        CREATEFIELD(edgeScalar,EdgeScalar,MESHTYPE,DEEP)                \
        CREATEFIELD(vertexScalar,VertexScalar,MESHTYPE,DEEP)            \
        CREATEFIELD(faceVector,FaceVector,MESHTYPE,DEEP)                \
        CREATEFIELD(edgeVector,EdgeVector,MESHTYPE,DEEP)                \
        CREATEFIELD(vertexVector,VertexVector,MESHTYPE,DEEP)            \
                                                                        \
        forAll(fvMsh, l)                                                \
        {                                                               \
            io.writeNow<MESHTYPE>(l);                                   \
        }                                                               \
                                                                        \
        RESETFIELD(Label,MESHTYPE)                                      \
        RESETFIELD(Scalar,MESHTYPE)                                     \
        RESETFIELD(Vector,MESHTYPE)                                     \
        RESETFIELD(Tensor,MESHTYPE)                                     \
        RESETFIELD(DiagTensor,MESHTYPE)                                 \
        RESETFIELD(SymmTensor,MESHTYPE)                                 \
        RESETFIELD(SphericalTensor,MESHTYPE)                            \
        RESETFIELD(FaceScalar,MESHTYPE)                                 \
        RESETFIELD(EdgeScalar,MESHTYPE)                                 \
        RESETFIELD(VertexScalar,MESHTYPE)                               \
        RESETFIELD(FaceVector,MESHTYPE)                                 \
        RESETFIELD(EdgeVector,MESHTYPE)                                 \
        RESETFIELD(VertexVector,MESHTYPE)                               \
                                                                        \
        forAll(fvMsh, l)                                                \
        {                                                               \
            io.read<MESHTYPE>(l);                                       \
        }                                                               \
                                                                        \
        TESTFIELD(label,Label,MESHTYPE,GHOSTS)                          \
        TESTFIELD(scalar,Scalar,MESHTYPE,GHOSTS)                        \
        TESTFIELD(vector,Vector,MESHTYPE,GHOSTS)                        \
        TESTFIELD(tensor,Tensor,MESHTYPE,GHOSTS)                        \
        TESTFIELD(diagTensor,DiagTensor,MESHTYPE,GHOSTS)                \
        TESTFIELD(symmTensor,SymmTensor,MESHTYPE,GHOSTS)                \
        TESTFIELD(sphericalTensor,SphericalTensor,MESHTYPE,GHOSTS)      \
        TESTFIELD(faceScalar,FaceScalar,MESHTYPE,GHOSTS)                \
        TESTFIELD(edgeScalar,EdgeScalar,MESHTYPE,GHOSTS)                \
        TESTFIELD(vertexScalar,VertexScalar,MESHTYPE,GHOSTS)            \
        TESTFIELD(faceVector,FaceVector,MESHTYPE,GHOSTS)                \
        TESTFIELD(edgeVector,EdgeVector,MESHTYPE,GHOSTS)                \
        TESTFIELD(vertexVector,VertexVector,MESHTYPE,GHOSTS)            \
    }

    // Test colocated

    TEST(colocated,1,false,false,true)
    TEST(colocated,2,true,false,true)
    TEST(colocated,3,false,true,true)
    TEST(colocated,4,true,true,true)
    TEST(colocated,1,false,false,false)
    TEST(colocated,2,true,false,false)
    TEST(colocated,3,false,true,false)
    TEST(colocated,4,true,true,false)

    // Test staggered only on structures mesh

    if (fvMsh.structured())
    {
        TEST(staggered,5,false,false,true)
        TEST(staggered,6,true,false,true)
        TEST(staggered,7,false,true,true)
        TEST(staggered,8,true,true,true)
        TEST(staggered,5,false,false,false)
        TEST(staggered,6,true,false,false)
        TEST(staggered,7,false,true,false)
        TEST(staggered,8,true,true,false)
    }
}
