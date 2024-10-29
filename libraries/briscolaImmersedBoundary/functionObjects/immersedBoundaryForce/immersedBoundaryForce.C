#include "immersedBoundaryForce.H"
#include "addToRunTimeSelectionTable.H"
#include "colocated.H"
#include "staggered.H"
#include "linearSystem.H"
#include "OSspecific.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

typedef immersedBoundaryForce<colocated> coloImmersedBoundaryForce;
typedef immersedBoundaryForce<staggered> stagImmersedBoundaryForce;

defineTemplateTypeNameAndDebugWithName
(
    coloImmersedBoundaryForce,
    "coloImmersedBoundaryForce",
    0
);

defineTemplateTypeNameAndDebugWithName
(
    stagImmersedBoundaryForce,
    "stagImmersedBoundaryForce",
    0
);

addToRunTimeSelectionTable
(
    functionObject,
    coloImmersedBoundaryForce,
    dictionary
);

addToRunTimeSelectionTable
(
    functionObject,
    stagImmersedBoundaryForce,
    dictionary
);

template<class MeshType>
immersedBoundaryForce<MeshType>::immersedBoundaryForce
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    ibName_(dict.lookup("ibName")),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    boundaryPtr_(nullptr)
{
    forAll(fvMsh_.immersedBoundaries<MeshType>(), i)
        if (fvMsh_.immersedBoundaries<MeshType>()[i].name() == ibName_)
            boundaryPtr_ = &fvMsh_.immersedBoundaries<MeshType>()[i];

    if (boundaryPtr_ == nullptr)
        FatalError
            << "Immersed boundary " << ibName_ << " requested but not found."
            << endl << abort(FatalError);

    if (Pstream::master())
    {
        const fileName path("postProcessing"/name_/runTime_.timeName());
        mkDir(path);

        OFstream file();

        filePtr_.reset(new OFstream(path/this->typeName + ".txt"));

        filePtr_()
            << "# time force" << endl;
    }
}

template<class MeshType>
immersedBoundaryForce<MeshType>::~immersedBoundaryForce()
{}

template<class MeshType>
bool immersedBoundaryForce<MeshType>::execute()
{
    typedef linearSystem<stencil,typename MeshType::vectorType,MeshType>
        sysType;

    typedef meshField<typename MeshType::vectorType,MeshType>
        fieldType;

    const sysType& sys =
        fvMsh_.db().lookupObject<sysType>
        (
            sysType::typeName + "(" + UName_ + ")"
        );

    fieldType r(- sys.template residual<true>());

    forAll(sys.x().immersedBoundaryConditions(), i)
    {
        const auto& ibc = sys.x().immersedBoundaryConditions()[i];

        if (&ibc.ib() == boundaryPtr_)
        {
            r *= ibc.forcingMask();
            break;
        }
    }

    force_ = gSum(r);

    if (Pstream::master())
    {
        if (log)
            Info<< "IBM force on " << ibName_ << " = " << force_ << nl << endl;

        filePtr_()
            << runTime_.time().value() << " " << force_ << endl;
    }

    return true;
}

template<class MeshType>
bool immersedBoundaryForce<MeshType>::write()
{
    return true;
}

template<class MeshType>
bool immersedBoundaryForce<MeshType>::end()
{
    return true;
}

}

}

}

}
