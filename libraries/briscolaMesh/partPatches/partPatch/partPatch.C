#include "partPatch.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(partPatch, 0);
defineRunTimeSelectionTable(partPatch, dictionary);

partPatch::partPatch(const mesh& msh, const dictionary& dict)
:
    msh_(msh),
    dict_(dict),
    name_(dict.lookup("name")),
    boundaryOffset_(dict.lookup("boundaryOffset")),
    T_(dict.lookupOrDefault<labelTensor>("T", eye)),
    master_(true)
{}

partPatch::partPatch
(
    const partPatch& pp
)
:
    msh_(pp.msh_),
    dict_(pp.dict_),
    name_(pp.name_),
    boundaryOffset_(pp.boundaryOffset_),
    T_(pp.T_),
    master_(pp.master_)
{}

partPatch::~partPatch()
{}

autoPtr<partPatch> partPatch::New
(
    const mesh& msh,
    const dictionary& dict
)
{
    const word partPatchType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(partPatchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown part patch type " << partPatchType << endl
            << "Valid partPatch types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<partPatch>(cstrIter()(msh, dict));
}

}

}
