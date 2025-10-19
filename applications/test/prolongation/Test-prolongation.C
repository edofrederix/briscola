#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "prolongationSchemes.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testConstantProlongation(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<prolongationScheme<Type,MeshType>> P
    (
        prolongationScheme<Type,MeshType>::NewType
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    f = Zero;

    // Check if a constant value is prolongated correctly

    for (label l = f.size()-1; l > 0; l--)
    {
        f[l] = pTraits<Type>::one*2.0;

        forAll(f[l], d)
        {
            P->prolong(f[l-1][d], f[l][d], eqOp<Type>());

            forAllCells(f[l-1][d], i, j, k)
                if (mag(f[l-1][d](i,j,k) - pTraits<Type>::one*2.0) > 1e-14)
                    FatalErrorInFunction
                        << "test 1 failed" << abort(FatalError);
        }
    }
}

template<class Type, class MeshType>
void testOneLinearProlongation(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<prolongationScheme<Type,MeshType>> P
    (
        prolongationScheme<Type,MeshType>::NewType
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

        for (label l = f.size()-1; l > 0; l--)
        {
            // Also set ghost values which are needed for prolongation

            forAll(f[l], d)
            for (label i = f.I(l,d).left()  -1; i < f.I(l,d).right()+1; i++)
            for (label j = f.I(l,d).bottom()-1; j < f.I(l,d).top()  +1; j++)
            for (label k = f.I(l,d).aft()   -1; k < f.I(l,d).fore() +1; k++)
            {
                f(l,d,i,j,k) = cc(l,d,i,j,k)[dir]*pTraits<Type>::one;
            }

            forAll(f[l], d)
            {
                P->prolong(f[l-1][d], f[l][d], eqOp<Type>());

                forAllCells(f[l-1][d], i, j, k)
                {
                    if
                    (
                        mag
                        (
                            f[l-1][d](i,j,k)
                          - cc[l-1][d](i,j,k)[dir]*pTraits<Type>::one
                        ) > 1e-14
                    )
                    {
                        FatalErrorInFunction
                            << "test 2 failed" << abort(FatalError);
                    }
                }
            }
        }
    }
}

template<class Type, class MeshType>
void testTwoLinearProlongation(const fvMesh& fvMsh, const word scheme)\
{
    autoPtr<prolongationScheme<Type,MeshType>> P
    (
        prolongationScheme<Type,MeshType>::NewType
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

        for (label l = f.size()-1; l > 0; l--)
        {
            // Also set ghost values which are needed for prolongation

            forAll(f[l], d)
            for (label i = f.I(l,d).left()  -1; i < f.I(l,d).right()+1; i++)
            for (label j = f.I(l,d).bottom()-1; j < f.I(l,d).top()  +1; j++)
            for (label k = f.I(l,d).aft()   -1; k < f.I(l,d).fore() +1; k++)
            {
                f(l,d,i,j,k) =
                    (cc(l,d,i,j,k)[dir1]+cc(l,d,i,j,k)[dir2])
                  * pTraits<Type>::one;
            }

            forAll(f[l], d)
            {
                P->prolong(f[l-1][d], f[l][d], eqOp<Type>());

                forAllCells(f[l-1][d], i, j, k)
                {
                    if
                    (
                        Pstream::myProcNo() == 0 &&
                        mag
                        (
                            f[l-1][d](i,j,k)
                          - (cc[l-1][d](i,j,k)[dir1] + cc[l-1][d](i,j,k)[dir2])
                          * pTraits<Type>::one
                        ) > 1e-14
                    )
                    {
                        FatalErrorInFunction
                            << "test 3 failed" << abort(FatalError);
                    }
                }
            }
        }
    }
}

template<class Type, class MeshType>
void testThreeLinearProlongation(const fvMesh& fvMsh, const word scheme)
{
    autoPtr<prolongationScheme<Type,MeshType>> P
    (
        prolongationScheme<Type,MeshType>::NewType
        (
            fvMsh,
            scheme
        )
    );

    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    f = Zero;

    for (label l = f.size()-1; l > 0; l--)
    {
        // Also set ghost values which are needed for prolongation

        forAll(f[l], d)
        for (label i = f.I(l,d).left()  -1; i < f.I(l,d).right()+1; i++)
        for (label j = f.I(l,d).bottom()-1; j < f.I(l,d).top()  +1; j++)
        for (label k = f.I(l,d).aft()   -1; k < f.I(l,d).fore() +1; k++)
        {
            f(l,d,i,j,k) =
                (
                    cc(l,d,i,j,k).x()
                  + cc(l,d,i,j,k).y()
                  + cc(l,d,i,j,k).z()
                )
              * pTraits<Type>::one;
        }

        forAll(f[l], d)
        {
            P->prolong(f[l-1][d], f[l][d], eqOp<Type>());

            forAllCells(f[l-1][d], i, j, k)
            {
                if
                (
                    Pstream::myProcNo() == 0 &&
                    mag
                    (
                        f[l-1][d](i,j,k)
                      - (
                            cc[l-1][d](i,j,k).x()
                          + cc[l-1][d](i,j,k).y()
                          + cc[l-1][d](i,j,k).z()
                        )
                      * pTraits<Type>::one
                    ) > 1e-14
                )
                {
                    FatalErrorInFunction
                        << "test 4 failed" << abort(FatalError);
                }
            }
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

    // Prolongation of a constant value

    testConstantProlongation<scalar,colocated>(fvMsh, "injection");
    testConstantProlongation<vector,colocated>(fvMsh, "injection");
    testConstantProlongation<tensor,colocated>(fvMsh, "injection");

    testConstantProlongation<scalar,staggered>(fvMsh, "injection");
    testConstantProlongation<vector,staggered>(fvMsh, "injection");
    testConstantProlongation<tensor,staggered>(fvMsh, "injection");

    testConstantProlongation<scalar,colocated>(fvMsh, "linear");
    testConstantProlongation<vector,colocated>(fvMsh, "linear");
    testConstantProlongation<tensor,colocated>(fvMsh, "linear");

    testConstantProlongation<scalar,staggered>(fvMsh, "linear");
    testConstantProlongation<vector,staggered>(fvMsh, "linear");
    testConstantProlongation<tensor,staggered>(fvMsh, "linear");

    // Prolongation of a linearly increasing value in one dimension. Does not
    // work for injection.

    testOneLinearProlongation<scalar,colocated>(fvMsh, "linear");
    testOneLinearProlongation<vector,colocated>(fvMsh, "linear");
    testOneLinearProlongation<tensor,colocated>(fvMsh, "linear");

    testOneLinearProlongation<scalar,staggered>(fvMsh, "linear");
    testOneLinearProlongation<vector,staggered>(fvMsh, "linear");
    testOneLinearProlongation<tensor,staggered>(fvMsh, "linear");

    // Prolongation of a linearly increasing value in two dimensions. Does not
    // work for injection.

    testTwoLinearProlongation<scalar,colocated>(fvMsh, "linear");
    testTwoLinearProlongation<vector,colocated>(fvMsh, "linear");
    testTwoLinearProlongation<tensor,colocated>(fvMsh, "linear");

    testTwoLinearProlongation<scalar,staggered>(fvMsh, "linear");
    testTwoLinearProlongation<vector,staggered>(fvMsh, "linear");
    testTwoLinearProlongation<tensor,staggered>(fvMsh, "linear");

    // Prolongation of a linearly increasing value in three dimensions. Does not
    // work for injection.

    testThreeLinearProlongation<scalar,colocated>(fvMsh, "linear");
    testThreeLinearProlongation<vector,colocated>(fvMsh, "linear");
    testThreeLinearProlongation<tensor,colocated>(fvMsh, "linear");

    testThreeLinearProlongation<scalar,staggered>(fvMsh, "linear");
    testThreeLinearProlongation<vector,staggered>(fvMsh, "linear");
    testThreeLinearProlongation<tensor,staggered>(fvMsh, "linear");
}
