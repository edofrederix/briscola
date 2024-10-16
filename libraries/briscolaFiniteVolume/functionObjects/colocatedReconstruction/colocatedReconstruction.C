#include "colocatedReconstruction.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

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
    reconstructedFields_(),
    requireBC_(dict.lookupOrDefault<Switch>("requireBC", false))
{
    init();
}

void colocatedReconstruction::init()
{
    const objectRegistry& db = fvMsh_.db();

    // Number of staggered scalar fields needed
    scalar size = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            size += 1;
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[i] << " requested for sampling by "
                << this->name() << " but not found in registry." << endl;
        }
    }

    reconstructedFields_.setSize(size);

    forAll(reconstructedFields_, i)
    {
        reconstructedFields_.set
        (
            i,
            new colocatedVectorField
            (
                fields_[i]+"c",
                fvMsh_,
                requireBC_ ? IOobject::MUST_READ : IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                true,
                requireBC_
            )
        );

        reconstructedFields_[i] = Zero;
    }
}

colocatedReconstruction::~colocatedReconstruction()
{}

bool colocatedReconstruction::execute()
{
    const objectRegistry& db = fvMsh_.db();

    label index = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            const staggeredScalarField& stagField
                = db.lookupObject<staggeredScalarField>(fields_[i]);

            reconstructedFields_[index] = ex::reconstruct(stagField);
            reconstructedFields_[index].correctBoundaryConditions();

            index++;
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
