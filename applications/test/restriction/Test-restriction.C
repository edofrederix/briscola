#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "restrictionSchemes.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testConstantRestriction(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<restrictionScheme<Type,MeshType>> S
    (
        restrictionScheme<Type,MeshType>::New
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    f = Zero;

    f[0] = pTraits<Type>::one*2.0;

    S->restrict(f);

    // Check if a constant value is restricted correctly

    forAllCells(f, l, d, i, j, k)
    if (l > 0)
        if (mag(f(l,d,i,j,k) - pTraits<Type>::one*2.0) > 1e-14)
            FatalErrorInFunction
                << "test 1 failed" << abort(FatalError);
}

template<class Type, class MeshType>
void testOneLinearRestriction(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<restrictionScheme<Type,MeshType>> S
    (
        restrictionScheme<Type,MeshType>::New
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    for (label dir = 0; dir < 3; dir++)
    {
        f = Zero;

        forAllCells(f, d, i, j, k)
        {
            f(d,i,j,k) = cc(d,i,j,k)[dir]*pTraits<Type>::one;
        }

        S->restrict(f);

        forAllCells(f, l, d, i, j, k)
        if (l > 0)
        {
            if
            (
                mag(f(l,d,i,j,k) - cc(l,d,i,j,k)[dir]*pTraits<Type>::one)
              > 1e-14
            )
            {
                FatalErrorInFunction
                    << "test 2 failed" << abort(FatalError);
            }
        }
    }
}

template<class Type, class MeshType>
void testTwoLinearRestriction(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<restrictionScheme<Type,MeshType>> S
    (
        restrictionScheme<Type,MeshType>::New
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    for (label dir1 = 0; dir1 < 3; dir1++)
    for (label dir2 = dir1+1; dir2 < 3; dir2++)
    {
        f = Zero;

        forAllCells(f, d, i, j, k)
        {
            f(d,i,j,k) =
                (cc(d,i,j,k)[dir1]+cc(d,i,j,k)[dir2])
              * pTraits<Type>::one;
        }

        S->restrict(f);

        forAllCells(f, l, d, i, j, k)
        if (l > 0)
        {
            if
            (
                mag
                (
                    f(l,d,i,j,k)
                  - (cc(l,d,i,j,k)[dir1] + cc(l,d,i,j,k)[dir2])
                  * pTraits<Type>::one
                )
              > 1e-14
            )
            {
                FatalErrorInFunction
                    << "test 3 failed" << abort(FatalError);
            }
        }
    }
}

template<class Type, class MeshType>
void testThreeLinearRestriction(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<restrictionScheme<Type,MeshType>> S
    (
        restrictionScheme<Type,MeshType>::New
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    f = Zero;

    forAllCells(f, d, i, j, k)
    {
        f(d,i,j,k) =
            (
                cc(d,i,j,k).x()
              + cc(d,i,j,k).y()
              + cc(d,i,j,k).z()
            )
          * pTraits<Type>::one;
    }

    S->restrict(f);

    forAllCells(f, l, d, i, j, k)
    if (l > 0)
    {
        if
        (
            mag
            (
                f(l,d,i,j,k)
              - (
                    cc(l,d,i,j,k).x()
                  + cc(l,d,i,j,k).y()
                  + cc(l,d,i,j,k).z()
                )
              * pTraits<Type>::one
            )
          > 1e-14
        )
        {
            FatalErrorInFunction
                << "test 4 failed" << abort(FatalError);
        }
    }
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

    // Restriction of a constant value

    testConstantRestriction<scalar,colocated>(fvMsh, "average");
    testConstantRestriction<vector,colocated>(fvMsh, "average");
    testConstantRestriction<tensor,colocated>(fvMsh, "average");

    testConstantRestriction<scalar,staggered>(fvMsh, "average");
    testConstantRestriction<vector,staggered>(fvMsh, "average");
    testConstantRestriction<tensor,staggered>(fvMsh, "average");

    testConstantRestriction<scalar,colocated>(fvMsh, "linear");
    testConstantRestriction<vector,colocated>(fvMsh, "linear");
    testConstantRestriction<tensor,colocated>(fvMsh, "linear");

    testConstantRestriction<scalar,staggered>(fvMsh, "linear");
    testConstantRestriction<vector,staggered>(fvMsh, "linear");
    testConstantRestriction<tensor,staggered>(fvMsh, "linear");

    // Restriction of a linearly increasing value in one dimension

    testOneLinearRestriction<scalar,colocated>(fvMsh, "average");
    testOneLinearRestriction<vector,colocated>(fvMsh, "average");
    testOneLinearRestriction<tensor,colocated>(fvMsh, "average");

    testOneLinearRestriction<scalar,staggered>(fvMsh, "average");
    testOneLinearRestriction<vector,staggered>(fvMsh, "average");
    testOneLinearRestriction<tensor,staggered>(fvMsh, "average");

    testOneLinearRestriction<scalar,colocated>(fvMsh, "linear");
    testOneLinearRestriction<vector,colocated>(fvMsh, "linear");
    testOneLinearRestriction<tensor,colocated>(fvMsh, "linear");

    testOneLinearRestriction<scalar,staggered>(fvMsh, "linear");
    testOneLinearRestriction<vector,staggered>(fvMsh, "linear");
    testOneLinearRestriction<tensor,staggered>(fvMsh, "linear");

    // Restriction of a linearly increasing value in two dimensions

    testTwoLinearRestriction<scalar,colocated>(fvMsh, "average");
    testTwoLinearRestriction<vector,colocated>(fvMsh, "average");
    testTwoLinearRestriction<tensor,colocated>(fvMsh, "average");

    testTwoLinearRestriction<scalar,staggered>(fvMsh, "average");
    testTwoLinearRestriction<vector,staggered>(fvMsh, "average");
    testTwoLinearRestriction<tensor,staggered>(fvMsh, "average");

    testTwoLinearRestriction<scalar,colocated>(fvMsh, "linear");
    testTwoLinearRestriction<vector,colocated>(fvMsh, "linear");
    testTwoLinearRestriction<tensor,colocated>(fvMsh, "linear");

    testTwoLinearRestriction<scalar,staggered>(fvMsh, "linear");
    testTwoLinearRestriction<vector,staggered>(fvMsh, "linear");
    testTwoLinearRestriction<tensor,staggered>(fvMsh, "linear");

    // Restriction of a linearly increasing value in three dimensions

    testThreeLinearRestriction<scalar,colocated>(fvMsh, "average");
    testThreeLinearRestriction<vector,colocated>(fvMsh, "average");
    testThreeLinearRestriction<tensor,colocated>(fvMsh, "average");

    testThreeLinearRestriction<scalar,staggered>(fvMsh, "average");
    testThreeLinearRestriction<vector,staggered>(fvMsh, "average");
    testThreeLinearRestriction<tensor,staggered>(fvMsh, "average");

    testThreeLinearRestriction<scalar,colocated>(fvMsh, "linear");
    testThreeLinearRestriction<vector,colocated>(fvMsh, "linear");
    testThreeLinearRestriction<tensor,colocated>(fvMsh, "linear");

    testThreeLinearRestriction<scalar,staggered>(fvMsh, "linear");
    testThreeLinearRestriction<vector,staggered>(fvMsh, "linear");
    testThreeLinearRestriction<tensor,staggered>(fvMsh, "linear");
}
