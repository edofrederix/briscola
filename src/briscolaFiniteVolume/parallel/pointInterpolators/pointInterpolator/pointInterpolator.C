#include "pointInterpolator.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(pointInterpolator<colocated>, 0);
defineTemplateTypeNameAndDebug(pointInterpolator<staggered>, 0);

defineTemplateRunTimeSelectionTable(pointInterpolator<colocated>, dictionary);
defineTemplateRunTimeSelectionTable(pointInterpolator<staggered>, dictionary);

template<class MeshType>
pointInterpolator<MeshType>::pointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points,
    const bool global,
    const label l,
    const label d
)
:
    fvMsh_(fvMsh),
    global_(global),
    start_(0),
    size_(points.size()),
    l_(l),
    d_(d)
{
    const meshDirection<vertexVector,MeshType>& v =
        fvMsh.template metrics<MeshType>().vertexCenters()[l][d];

    if (!global)
    {
        List<vectorList> globalPoints(Pstream::nProcs());
        globalPoints[Pstream::myProcNo()] = points;

        Pstream::gatherList
        (
            globalPoints,
            Pstream::msgType(),
            fvMsh[l].comms()
        );

        Pstream::scatterList
        (
            globalPoints,
            Pstream::msgType(),
            fvMsh[l].comms()
        );

        int nPoints = 0;
        forAll(globalPoints, i)
            nPoints += globalPoints[i].size();

        points_.setSize(nPoints);

        int start = 0;
        forAll(globalPoints, i)
        {
            SubList<vector> sub(points_, globalPoints[i].size(), start);
            sub = globalPoints[i];

            if (i == Pstream::myProcNo())
                start_ = start;

            start += globalPoints[i].size();
        }
    }
    else
    {
        points_ = points;
    }

    indices_.setSize(points_.size());
    cellCoordinates_.setSize(points_.size());

    forAll(points_, i)
    {
        indices_[i] = fvMsh.findCell<MeshType>(trimPrecision(points_[i]), l, d);

        if (indices_[i] != -unitXYZ)
        {
            cellCoordinates_[i] =
                interpolationWeights(points_[i], v(indices_[i]));
        }
        else
        {
            cellCoordinates_[i] = -vector::one;
        }
    }
}

template<class MeshType>
pointInterpolator<MeshType>::pointInterpolator
(
    const pointInterpolator<MeshType>& interp
)
:
    fvMsh_(interp.fvMsh_),
    global_(interp.global_),
    l_(interp.l_),
    d_(interp.d_),
    points_(interp.points_),
    indices_(interp.indices_),
    cellCoordinates_(interp.cellCoordinates_)
{}

template<class MeshType>
pointInterpolator<MeshType>::~pointInterpolator()
{}

template<class MeshType>
autoPtr<pointInterpolator<MeshType>> pointInterpolator<MeshType>::New
(
    const fvMesh& fvMsh,
    const vectorList& points,
    const word type,
    const bool global,
    const label l,
    const label d
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown point interpolator type " << type
            << ". Valid point interpolators types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pointInterpolator<MeshType>>
    (
        cstrIter()(fvMsh,points,global,l,d)
    );
}

template<class MeshType>
vectorList pointInterpolator<MeshType>::missingPoints() const
{
    List<label> found(indices_.size());

    forAll(indices_, i)
        found[i] = (indices_[i] != -unitXYZ);

    Pstream::listCombineGather
    (
        found,
        plusEqOp<label>(),
        Pstream::msgType(),
        fvMsh_[l_].comms()
    );

    DynamicList<vector> missing;

    if (Pstream::master())
    {
        forAll(found, i)
            if (!found[i])
                missing.append(points_[i]);

        for (int proc = 0; proc < Pstream::nProcs(); proc++)
        if (proc != Pstream::masterNo())
        {
            OPstream send
            (
                Pstream::commsTypes::blocking,
                proc,
                0,
                UPstream::msgType(),
                fvMsh_[l_].comms()
            );

            send << missing;
        }
    }
    else
    {
        IPstream recv
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo(),
            0,
            UPstream::msgType(),
            fvMsh_[l_].comms()
        );

        recv >> missing;
    }

    return missing;
}

template class pointInterpolator<colocated>;
template class pointInterpolator<staggered>;

}

}

}
