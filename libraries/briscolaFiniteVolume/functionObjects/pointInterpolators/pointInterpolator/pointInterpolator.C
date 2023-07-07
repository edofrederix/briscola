#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(pointInterpolator, 0);
defineRunTimeSelectionTable(pointInterpolator, dictionary);

pointInterpolator::pointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points
)
:
    fvMsh_(fvMsh),
    points_(points),
    indices_(points.size()),
    cellCoordinates_(points.size())
{
    const mesh& msh = fvMsh.msh();

    forAll(points, i)
    {
        indices_[i] = msh.findCell(points[i]);

        if (indices_[i] != -unitXYZ)
        {
            cellCoordinates_[i] =
                msh[0].points().cellCoordinates(points[i], indices_[i], true);
        }
        else
        {
            cellCoordinates_ = -vector::one;
        }
    }
}

pointInterpolator::pointInterpolator(const pointInterpolator& interp)
:
    fvMsh_(interp.fvMsh_),
    points_(interp.points_),
    indices_(interp.indices_),
    cellCoordinates_(interp.cellCoordinates_)
{}

pointInterpolator::~pointInterpolator()
{}

autoPtr<pointInterpolator> pointInterpolator::New
(
    const fvMesh& fvMsh,
    const vectorList& points,
    const word type
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown point interpolator type " << type
            << ". Valid point interpolators types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pointInterpolator>(cstrIter()(fvMsh, points));
}

bool pointInterpolator::allPointsFound() const
{
    List<label> found(indices_.size());

    forAll(indices_, i)
        found[i] = (indices_[i] != -unitXYZ);

    Pstream::listCombineGather(found, plusEqOp<label>());
    Pstream::listCombineScatter(found);

    forAll(found, i)
        if (!found[i])
            return false;

    return true;
}

}

}

}
