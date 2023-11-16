#include "partPatch.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(partPatch, 0);
defineRunTimeSelectionTable(partPatch, dictionary);

void partPatch::extend()
{
    for (int i = 0; i < 6; i++)
        if (boundaryOffset_[i/2] == 0)
            extension_[i] = 1;
}

partPatch::partPatch(const mesh& msh, const dictionary& dict)
:
    dict_(dict),
    name_(dict.lookup("name")),
    boundaryOffset_(dict.lookup("boundaryOffset")),
    boundaryOffsetDegree_
    (
        cmptSum(cmptMag(boundaryOffset_))
    ),
    T_(dict.lookupOrDefault<labelTensor>("T", eye)),
    master_(true),
    extension_(Zero)
{}

partPatch::partPatch
(
    const partPatch& pp
)
:
    dict_(pp.dict_),
    name_(pp.name_),
    boundaryOffset_(pp.boundaryOffset_),
    boundaryOffsetDegree_(pp.boundaryOffsetDegree_),
    T_(pp.T_),
    master_(pp.master_),
    extension_(Zero)
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
