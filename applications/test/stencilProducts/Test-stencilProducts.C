#include "stencil.H"
#include "symmStencil.H"
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
    meshField<SType,MeshType> f("f", fvMsh);
    meshField<Type,MeshType> x("x", fvMsh);

    f = pTraits<SType>::one;
    x = pTraits<Type>::one;

    Type s(pTraits<Type>::one);

    // Fields

    Amul(f,x);
    Amul(f,(x*2.0));
    Amul((2.0*f),x);
    Amul((2.0*f),(x*2.0));

    Amul(x,f);
    Amul((x*2.0),f);
    Amul(x,(2.0*f));
    Amul((x*2.0),(2.0*f));

    Amul(f,s);
    Amul(s,f);
    Amul((2.0*f),s);
    Amul(s,(2.0*f));

    rowSum(f);
    rowSum(2.0*f);

    matrixSum(f);
    matrixSum(2.0*f);

    neighborSum(f);
    neighborSum(2.0*f);

    // Levels

    Amul(f[0],x[0]);
    Amul(f[0],(x[0]*2.0));
    Amul((2.0*f[0]),x[0]);
    Amul((2.0*f[0]),(x[0]*2.0));

    Amul(x[0],f[0]);
    Amul((x[0]*2.0),f[0]);
    Amul(x[0],(2.0*f[0]));
    Amul((x[0]*2.0),(2.0*f[0]));

    Amul(f[0],s);
    Amul(s,f[0]);
    Amul((2.0*f[0]),s);
    Amul(s,(2.0*f[0]));

    rowSum(f[0]);
    rowSum(2.0*f[0]);

    matrixSum(f[0]);
    matrixSum(2.0*f[0]);

    neighborSum(f[0]);
    neighborSum(2.0*f[0]);

    // Directions

    Amul(f.direction(),x.direction());
    Amul(f.direction(),(x.direction()*2.0));
    Amul((2.0*f.direction()),x.direction());
    Amul((2.0*f.direction()),(x.direction()*2.0));

    Amul(x.direction(),f.direction());
    Amul((x.direction()*2.0),f.direction());
    Amul(x.direction(),(2.0*f.direction()));
    Amul((x.direction()*2.0),(2.0*f.direction()));

    Amul(f.direction(),s);
    Amul(s,f.direction());
    Amul((2.0*f.direction()),s);
    Amul(s,(2.0*f.direction()));

    rowSum(f.direction());
    rowSum(2.0*f.direction());

    matrixSum(f.direction());
    matrixSum(2.0*f.direction());

    neighborSum(f.direction());
    neighborSum(2.0*f.direction());
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
