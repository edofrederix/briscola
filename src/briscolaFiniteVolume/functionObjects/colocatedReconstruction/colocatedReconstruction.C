#include "colocatedReconstruction.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemesReconstruction.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(colocatedReconstruction, 0);

addToRunTimeSelectionTable
(
    functionObject,
    colocatedReconstruction,
    dictionary
);

colocatedReconstruction::colocatedReconstruction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    fields_(dict.lookup("fields")),
    vectorFields_(),
    tensorFields_()
{
    init();
}

void colocatedReconstruction::init()
{
    const objectRegistry& db = fvMsh_.db();

    forAll(fields_, i)
    {
        if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            vectorFields_.append
            (
                new colocatedVectorField
                (
                    fields_[i] + "c",
                    fvMsh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    true
                )
            );
        }
        else if (db.foundObject<staggeredVectorField>(fields_[i]))
        {
            tensorFields_.append
            (
                new colocatedTensorField
                (
                    fields_[i] + "c",
                    fvMsh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    true
                )
            );
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[i] << " requested for reconstruction by "
                << this->name() << " but not found in registry." << endl;
        }
    }

    // Initialize to zero

    forAll(vectorFields_, i)
        vectorFields_[i] = Zero;

    forAll(tensorFields_, i)
        tensorFields_[i] = Zero;
}

colocatedReconstruction::~colocatedReconstruction()
{}

bool colocatedReconstruction::execute()
{
    const objectRegistry& db = fvMsh_.db();

    label iv = 0;
    label it = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            const staggeredScalarField& field =
                db.lookupObject<staggeredScalarField>(fields_[i]);

            vectorFields_[iv] = ex::reconstruct(field);
            vectorFields_[iv].correctBoundaryConditions();

            iv++;
        }
        else if (db.foundObject<staggeredVectorField>(fields_[i]))
        {
            const staggeredVectorField& field =
                db.lookupObject<staggeredVectorField>(fields_[i]);

            tensorFields_[it] = ex::reconstruct(field);
            tensorFields_[it].correctBoundaryConditions();

            it++;
        }
    }

    return true;
}

bool colocatedReconstruction::write()
{
    return true;
}

bool colocatedReconstruction::end()
{
    return true;
}

}

}

}

}
