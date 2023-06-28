#include "stencil.H"
#include "diagStencil.H"

#include "arguments.H"
#include "Time.H"

#include "colocated.H"
#include "staggered.H"
#include "meshField.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
void testStencilProducts(const fvMesh& fvMsh)
{
    meshField<SType,MeshType> stencil("stencil", fvMsh);
    meshField<Type,MeshType> x("x", fvMsh);

    stencil = pTraits<SType>::one;
    x = pTraits<Type>::one;

    Type s(pTraits<Type>::one);

    // Fields

    Amul(stencil,x);
    Amul(stencil,(x*2.0));
    Amul((2.0*stencil),x);
    Amul((2.0*stencil),(x*2.0));

    Amul(x,stencil);
    Amul((x*2.0),stencil);
    Amul(x,(2.0*stencil));
    Amul((x*2.0),(2.0*stencil));

    Amul(stencil,s);
    Amul(s,stencil);
    Amul((2.0*stencil),s);
    Amul(s,(2.0*stencil));

    rowSum(stencil);
    rowSum(2.0*stencil);

    matrixSum(stencil);
    matrixSum(2.0*stencil);

    neighborSum(stencil);
    neighborSum(2.0*stencil);

    // Levels

    Amul(stencil[0],x[0]);
    Amul(stencil[0],(x[0]*2.0));
    Amul((2.0*stencil[0]),x[0]);
    Amul((2.0*stencil[0]),(x[0]*2.0));

    Amul(x[0],stencil[0]);
    Amul((x[0]*2.0),stencil[0]);
    Amul(x[0],(2.0*stencil[0]));
    Amul((x[0]*2.0),(2.0*stencil[0]));

    Amul(stencil[0],s);
    Amul(s,stencil[0]);
    Amul((2.0*stencil[0]),s);
    Amul(s,(2.0*stencil[0]));

    rowSum(stencil[0]);
    rowSum(2.0*stencil[0]);

    matrixSum(stencil[0]);
    matrixSum(2.0*stencil[0]);

    neighborSum(stencil[0]);
    neighborSum(2.0*stencil[0]);

    // Directions

    Amul(stencil[0][0],x[0][0]);
    Amul(stencil[0][0],(x[0][0]*2.0));
    Amul((2.0*stencil[0][0]),x[0][0]);
    Amul((2.0*stencil[0][0]),(x[0][0]*2.0));

    Amul(x[0][0],stencil[0][0]);
    Amul((x[0][0]*2.0),stencil[0][0]);
    Amul(x[0][0],(2.0*stencil[0][0]));
    Amul((x[0][0]*2.0),(2.0*stencil[0][0]));

    Amul(stencil[0][0],s);
    Amul(s,stencil[0][0]);
    Amul((2.0*stencil[0][0]),s);
    Amul(s,(2.0*stencil[0][0]));

    rowSum(stencil[0][0]);
    rowSum(2.0*stencil[0][0]);

    matrixSum(stencil[0][0]);
    matrixSum(2.0*stencil[0][0]);

    neighborSum(stencil[0][0]);
    neighborSum(2.0*stencil[0][0]);
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

    testStencilProducts<stencil,label,colocated>(fvMsh);
    testStencilProducts<stencil,scalar,colocated>(fvMsh);
    testStencilProducts<stencil,vector,colocated>(fvMsh);
    testStencilProducts<stencil,tensor,colocated>(fvMsh);
    testStencilProducts<stencil,symmTensor,colocated>(fvMsh);
    testStencilProducts<stencil,sphericalTensor,colocated>(fvMsh);
    testStencilProducts<stencil,diagTensor,colocated>(fvMsh);

    testStencilProducts<diagStencil,label,colocated>(fvMsh);
    testStencilProducts<diagStencil,scalar,colocated>(fvMsh);
    testStencilProducts<diagStencil,vector,colocated>(fvMsh);
    testStencilProducts<diagStencil,tensor,colocated>(fvMsh);
    testStencilProducts<diagStencil,symmTensor,colocated>(fvMsh);
    testStencilProducts<diagStencil,sphericalTensor,colocated>(fvMsh);
    testStencilProducts<diagStencil,diagTensor,colocated>(fvMsh);

    // Staggered

    testStencilProducts<stencil,label,staggered>(fvMsh);
    testStencilProducts<stencil,scalar,staggered>(fvMsh);
    testStencilProducts<stencil,vector,staggered>(fvMsh);
    testStencilProducts<stencil,tensor,staggered>(fvMsh);
    testStencilProducts<stencil,symmTensor,staggered>(fvMsh);
    testStencilProducts<stencil,sphericalTensor,staggered>(fvMsh);
    testStencilProducts<stencil,diagTensor,staggered>(fvMsh);

    testStencilProducts<diagStencil,label,staggered>(fvMsh);
    testStencilProducts<diagStencil,scalar,staggered>(fvMsh);
    testStencilProducts<diagStencil,vector,staggered>(fvMsh);
    testStencilProducts<diagStencil,tensor,staggered>(fvMsh);
    testStencilProducts<diagStencil,symmTensor,staggered>(fvMsh);
    testStencilProducts<diagStencil,sphericalTensor,staggered>(fvMsh);
    testStencilProducts<diagStencil,diagTensor,staggered>(fvMsh);
}
