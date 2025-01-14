#include "grading.H"
#include "brick.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(grading, 0);
defineRunTimeSelectionTable(grading, dictionary);

grading::grading(const brick& b)
:
    b_(b)
{}

grading::grading(const grading& g)
:
    b_(g.b_)
{}

grading::~grading()
{}

autoPtr<grading> grading::New(const brick& b)
{
    word gradingType("none");

    if (b.dict().found("grading"))
    {
        b.dict().lookup("grading") >> gradingType;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(gradingType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown grading type " << gradingType << " for " << b << endl
            << "Valid grading types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<grading>(cstrIter()(b));
}

}

}
