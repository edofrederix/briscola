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
    meshField<Type,MeshType> field
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ
    );

    forAllCells(field, l, d, i, j, k)
    {
        field(l,d,i,j,k) = scalar(i+j+k)*pTraits<Type>::one;
    }

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

        const labelVector bo(bcd.offset());

        if (bcd.offsetDegree() != 1)
        {
            continue;
        }

        forAll(field, l)
        forAll(field[l], d)
        {
            const labelVector S(field.fvMsh().template S<MeshType>(l,d,bo));
            const labelVector E(field.fvMsh().template E<MeshType>(l,d,bo));

            if (MeshType::shifted(d,bo))
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (Foam::mag(field(l,d,ijk+bo)) > 1e-12)
                        FatalErrorInFunction
                            << "Test 1a failed" << endl << abort(FatalError);
                }
            }
            else
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (Foam::mag(field(l,d,ijk+bo) + field(l,d,ijk)) > 1e-12)
                        FatalErrorInFunction
                            << "Test 1b failed" << endl << abort(FatalError);
                }
            }
        }
    }

    // Test global boundary condition base types

    for (int i = 0; i < 6; i++)
    {
        const labelVector bo = faceOffsets[i];

        if(globalBoundaryConditionBaseType(field, bo) != DIRICHLETBC)
        {
            FatalErrorInFunction
                << "Test 1b failed" << endl
                << bo << endl
                << abort(FatalError);
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

    forAllCells(field, l, d, i, j, k)
        field(l,d,i,j,k) = scalar(i+j+k)*pTraits<Type>::one;

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

        const labelVector bo(bcd.offset());

        if (bcd.offsetDegree() != 1)
        {
            continue;
        }

        forAll(field, l)
        forAll(field[l], d)
        {
            const labelVector S(field.fvMsh().template S<MeshType>(l,d,bo));
            const labelVector E(field.fvMsh().template E<MeshType>(l,d,bo));

            if (MeshType::shifted(d,bo))
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (field(l,d,ijk+bo) != field(l,d,ijk-bo))
                        FatalErrorInFunction
                            << "Test 2a failed" << endl << abort(FatalError);
                }
            }
            else
            {
                for (label i = S.x(); i < E.x(); i++)
                for (label j = S.y(); j < E.y(); j++)
                for (label k = S.z(); k < E.z(); k++)
                {
                    labelVector ijk(i,j,k);

                    if (field(l,d,ijk+bo) != field(l,d,ijk))
                        FatalErrorInFunction
                            << "Test 2b failed" << endl << abort(FatalError);
                }
            }
        }
    }

    // Test global boundary condition base types

    for (int i = 0; i < 6; i++)
    {
        const labelVector bo = faceOffsets[i];

        if(globalBoundaryConditionBaseType(field, bo) != NEUMANNBC)
        {
            FatalErrorInFunction
                << "Test 2c failed" << endl << abort(FatalError);
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

    testDirichlet<scalar,staggered>(fvMsh);
    testDirichlet<vector,staggered>(fvMsh);


    testNeumann<scalar,colocated>(fvMsh);
    testNeumann<vector,colocated>(fvMsh);
    testNeumann<tensor,colocated>(fvMsh);

    testNeumann<scalar,staggered>(fvMsh);
    testNeumann<vector,staggered>(fvMsh);
}
