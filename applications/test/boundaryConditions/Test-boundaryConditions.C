#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testDirichlet(const fvMesh& fvMsh)
{
    Type value = 2.0*pTraits<Type>::one;

    meshField<Type,MeshType> field
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ
    );

    forAll(field, l)
    forAll(field[l], d)
    forAllCells(field[l][d], i, j, k)
    {
        field[l][d](i,j,k) = scalar(i+j+k)*pTraits<Type>::one;
    }

    // We need to correct twice because of edges/vertices

    field.correctBoundaryConditions();
    field.correctBoundaryConditions();

    forAll(field.boundaryConditions(), bci)
    if
    (
        field.boundaryConditions()[bci].baseType()
     == boundaryConditionBaseType::DIRICHLETBC
    )
    {
        const boundaryCondition<Type,MeshType>& bcd =
            field.boundaryConditions()[bci];

        const labelVector bo(bcd.boundaryOffset());

        if (bcd.boundaryOffsetDegree() > 1)
        {
            continue;
        }

        forAll(field, l)
        forAll(field[l], d)
        {
            const meshDirection<Type,MeshType>& dd = field[l][d];

            const labelVector S(dd.boundaryStart(bo));
            const labelVector E(dd.boundaryEnd(bo));

            if (dd.shifted(bo))
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (dd(ijk+bo) != value)
                    {
                        FatalErrorInFunction
                            << "Test 1a failed" << endl << abort(FatalError);
                    }
                }
            }
            else
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (Pstream::master())
                    if (dd(ijk+bo) != 2.0*value - dd(ijk))
                    {
                        FatalErrorInFunction
                            << "Test 1b failed" << endl << abort(FatalError);
                    }
                }
            }
        }
    }
}

template<class Type, class MeshType>
void testNeumann(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> field
    (
        "g-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ
    );

    forAll(field, l)
    forAll(field[l], d)
    forAllCells(field[l][d], i, j, k)
    {
        field[l][d](i,j,k) = scalar(i+j+k)*pTraits<Type>::one;
    }

    // We need to correct twice because of edges/vertices

    field.correctBoundaryConditions();
    field.correctBoundaryConditions();

    forAll(field.boundaryConditions(), bci)
    if
    (
        field.boundaryConditions()[bci].baseType()
     == boundaryConditionBaseType::NEUMANNBC
    )
    {
        const boundaryCondition<Type,MeshType>& bcd =
            field.boundaryConditions()[bci];

        const labelVector bo(bcd.boundaryOffset());

        if (bcd.boundaryOffsetDegree() > 1)
        {
            continue;
        }

        forAll(field, l)
        forAll(field[l], d)
        {
            const meshDirection<Type,MeshType>& dd = field[l][d];

            const labelVector S(dd.boundaryStart(bo));
            const labelVector E(dd.boundaryEnd(bo));

            if (dd.shifted(bo))
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (dd(ijk+bo) != dd(ijk-bo))
                    {
                        FatalErrorInFunction
                            << "Test 2a failed" << endl << abort(FatalError);
                    }
                }
            }
            else
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (Pstream::master())
                    if (dd(ijk+bo) != dd(ijk))
                    {
                        FatalErrorInFunction
                            << "Test 2b failed" << endl << abort(FatalError);
                    }
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


    testDirichlet<scalar,colocated>(fvMsh);
    testDirichlet<vector,colocated>(fvMsh);

    testDirichlet<tensor,colocated>(fvMsh);
    testDirichlet<diagTensor,colocated>(fvMsh);
    testDirichlet<symmTensor,colocated>(fvMsh);
    testDirichlet<sphericalTensor,colocated>(fvMsh);

    testDirichlet<stencil,colocated>(fvMsh);
    testDirichlet<diagStencil,colocated>(fvMsh);

    testDirichlet<faceScalar,colocated>(fvMsh);
    testDirichlet<faceVector,colocated>(fvMsh);


    testDirichlet<scalar,staggered>(fvMsh);
    testDirichlet<vector,staggered>(fvMsh);

    testDirichlet<tensor,staggered>(fvMsh);
    testDirichlet<diagTensor,staggered>(fvMsh);
    testDirichlet<symmTensor,staggered>(fvMsh);
    testDirichlet<sphericalTensor,staggered>(fvMsh);

    testDirichlet<stencil,staggered>(fvMsh);
    testDirichlet<diagStencil,staggered>(fvMsh);

    testDirichlet<faceScalar,staggered>(fvMsh);
    testDirichlet<faceVector,staggered>(fvMsh);


    testNeumann<scalar,colocated>(fvMsh);
    testNeumann<vector,colocated>(fvMsh);

    testNeumann<tensor,colocated>(fvMsh);
    testNeumann<diagTensor,colocated>(fvMsh);
    testNeumann<symmTensor,colocated>(fvMsh);
    testNeumann<sphericalTensor,colocated>(fvMsh);

    testNeumann<stencil,colocated>(fvMsh);
    testNeumann<diagStencil,colocated>(fvMsh);

    testNeumann<faceScalar,colocated>(fvMsh);
    testNeumann<faceVector,colocated>(fvMsh);


    testNeumann<scalar,staggered>(fvMsh);
    testNeumann<vector,staggered>(fvMsh);

    testNeumann<tensor,staggered>(fvMsh);
    testNeumann<diagTensor,staggered>(fvMsh);
    testNeumann<symmTensor,staggered>(fvMsh);
    testNeumann<sphericalTensor,staggered>(fvMsh);

    testNeumann<stencil,staggered>(fvMsh);
    testNeumann<diagStencil,staggered>(fvMsh);

    testNeumann<faceScalar,staggered>(fvMsh);
    testNeumann<faceVector,staggered>(fvMsh);
}
