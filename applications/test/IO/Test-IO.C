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
    const fvMesh& fvMsh
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
            true
        )
    );

    meshField<Type,MeshType>& f = tf.ref();

    forAll(f, l)
    forAll(f[l], d)
    forAllCells(f[l][d], i, j, k)
    {
        f[l][d](i,j,k) = (Pstream::myProcNo()+i+j+k+l+d)*pTraits<Type>::one;
    }

    return tf;
}

template<class Type, class MeshType>
void testField(const meshField<Type,MeshType>& f)
{
    forAll(f, l)
    forAll(f[l], d)
    forAllCells(f[l][d], i, j, k)
    {
        if
        (
            f[l][d](i,j,k)
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
                << "expected value = "
                << (Pstream::myProcNo()+i+j+k)*pTraits<Type>::one << nl
                << "actual value = " << f[l][d](i,j,k)
                << abort(FatalError);
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaIO.H"

    runTime.setTime(1.0, 1);

    // Macros

    #define CREATEFIELD(TYPE, TYPENAME, MESHTYPE)                       \
        tmp<meshField<TYPE,MESHTYPE>> MESHTYPE##TYPENAME##TestField     \
        (                                                               \
            createField<TYPE,MESHTYPE>                                  \
            (                                                           \
                #MESHTYPE #TYPENAME "TestField",                        \
                fvMsh                                                   \
            )                                                           \
        );

    #define RESETFIELD(TYPENAME, MESHTYPE)                              \
        MESHTYPE##TYPENAME##TestField.ref() *= 0;                       \
        MESHTYPE##TYPENAME##TestField.ref().readOpt() =                 \
            IOobject::MUST_READ;

    #define TESTFIELD(TYPE, TYPENAME, MESHTYPE)                         \
        testField<TYPE,MESHTYPE>(MESHTYPE##TYPENAME##TestField.ref());

    // Colocated fields

    CREATEFIELD(label,Label,colocated)
    CREATEFIELD(scalar,Scalar,colocated)
    CREATEFIELD(vector,Vector,colocated)
    CREATEFIELD(tensor,Tensor,colocated)
    CREATEFIELD(diagTensor,DiagTensor,colocated)
    CREATEFIELD(symmTensor,SymmTensor,colocated)
    CREATEFIELD(sphericalTensor,SphericalTensor,colocated)
    CREATEFIELD(faceScalar,FaceScalar,colocated)
    CREATEFIELD(faceVector,FaceVector,colocated)

    forAll(fvMsh, l)
    {
        io.writeNow<colocated>(l);
    }

    RESETFIELD(Label,colocated)
    RESETFIELD(Scalar,colocated)
    RESETFIELD(Vector,colocated)
    RESETFIELD(Tensor,colocated)
    RESETFIELD(DiagTensor,colocated)
    RESETFIELD(SymmTensor,colocated)
    RESETFIELD(SphericalTensor,colocated)
    RESETFIELD(FaceScalar,colocated)
    RESETFIELD(FaceVector,colocated)

    forAll(fvMsh, l)
    {
        io.read<colocated>(l);
    }

    TESTFIELD(label,Label,colocated)
    TESTFIELD(scalar,Scalar,colocated)
    TESTFIELD(vector,Vector,colocated)
    TESTFIELD(tensor,Tensor,colocated)
    TESTFIELD(diagTensor,DiagTensor,colocated)
    TESTFIELD(symmTensor,SymmTensor,colocated)
    TESTFIELD(sphericalTensor,SphericalTensor,colocated)
    TESTFIELD(faceScalar,FaceScalar,colocated)
    TESTFIELD(faceVector,FaceVector,colocated)

    // Staggered fields only on structured meshes

    if (fvMsh.structured())
    {
        CREATEFIELD(label,Label,staggered)
        CREATEFIELD(scalar,Scalar,staggered)
        CREATEFIELD(vector,Vector,staggered)
        CREATEFIELD(tensor,Tensor,staggered)
        CREATEFIELD(diagTensor,DiagTensor,staggered)
        CREATEFIELD(symmTensor,SymmTensor,staggered)
        CREATEFIELD(sphericalTensor,SphericalTensor,staggered)
        CREATEFIELD(faceScalar,FaceScalar,staggered)
        CREATEFIELD(faceVector,FaceVector,staggered)

        forAll(fvMsh, l)
        {
            io.writeNow<staggered>(l);
        }

        RESETFIELD(Label,staggered)
        RESETFIELD(Scalar,staggered)
        RESETFIELD(Vector,staggered)
        RESETFIELD(Tensor,staggered)
        RESETFIELD(DiagTensor,staggered)
        RESETFIELD(SymmTensor,staggered)
        RESETFIELD(SphericalTensor,staggered)
        RESETFIELD(FaceScalar,staggered)
        RESETFIELD(FaceVector,staggered)

        forAll(fvMsh, l)
        {
            io.read<staggered>(l);
        }

        TESTFIELD(label,Label,staggered)
        TESTFIELD(scalar,Scalar,staggered)
        TESTFIELD(vector,Vector,staggered)
        TESTFIELD(tensor,Tensor,staggered)
        TESTFIELD(diagTensor,DiagTensor,staggered)
        TESTFIELD(symmTensor,SymmTensor,staggered)
        TESTFIELD(sphericalTensor,SphericalTensor,staggered)
        TESTFIELD(faceScalar,FaceScalar,staggered)
        TESTFIELD(faceVector,FaceVector,staggered)
    }
}
