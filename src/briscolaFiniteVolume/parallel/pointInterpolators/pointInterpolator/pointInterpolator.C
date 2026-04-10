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
    size_(points.size()),
    l_(l),
    d_(d),
    points_(points)
{
    const meshDirection<vertexVector,MeshType>& v =
        fvMsh.template metrics<MeshType>().vertexCenters()[l][d];

    // Find local points

    indices_.setSize(size_);
    cellCoordinates_.setSize(size_);
    requesters_.setSize(size_);
    providers_.setSize(size_);

    cellCoordinates_ = -vector::one;
    requesters_ = Pstream::myProcNo();
    providers_ = -1;

    forAll(points_, i)
    {
        indices_[i] =
            fvMsh.findCell<MeshType>(trimPrecision(points_[i]), l, d);

        if (indices_[i] != -unitXYZ)
        {
            cellCoordinates_[i] =
                interpolationWeights(points_[i], v(indices_[i]));

            providers_[i] = Pstream::myProcNo();
        }
    }

    if (global)
    {
        // If global then all processors have all points. Combine the providers
        // list to determine which processor provides which point.

        Pstream::listCombineGather(providers_, maxEqOp<label>());
        Pstream::listCombineScatter(providers_);
    }
    else
    {
        // For non-local points we need to find the corresponding processors

        DynamicList<vector> others;

        forAll(points_, i)
            if (providers_[i] == -1)
                others.append(points_[i]);

        List<vectorList> all(Pstream::nProcs());
        all[Pstream::myProcNo()] = others;

        // Distribute

        Pstream::gatherList(all);
        Pstream::scatterList(all);

        // Check others' points and when found set ourself as provider

        List<labelList> providers(Pstream::nProcs());
        forAll(providers, proc)
            providers[proc].resize(all[proc].size());

        forAll(all, proc)
        {
            providers[proc] = -1;

            forAll(all[proc], i)
            {
                const vector point = all[proc][i];

                const labelVector index =
                    fvMsh.findCell<MeshType>(trimPrecision(point), l, d);

                if (index != -unitXYZ)
                {
                    // If the point is ours, append it to the list of points to
                    // be interpolated

                    providers[proc][i] = Pstream::myProcNo();

                    points_.append(point);

                    requesters_.append(proc);
                    providers_.append(Pstream::myProcNo());

                    indices_.append(index);

                    cellCoordinates_.append
                    (
                        interpolationWeights(point, v(index))
                    );
                }
            }
        }

        // Distribute providers so that everyone is aware of who's doing what

        forAll(providers, proc)
        {
            Pstream::listCombineGather(providers[proc], maxEqOp<label>());
            Pstream::listCombineScatter(providers[proc]);
        }

        // Set providers

        label c = 0;
        forAll(providers_, i)
            if (providers_[i] == -1)
                providers_[i] = providers[Pstream::myProcNo()][c++];
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
    size_(interp.size_),
    l_(interp.l_),
    d_(interp.d_),
    points_(interp.points_),
    requesters_(interp.requesters_),
    providers_(interp.providers_),
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
    DynamicList<vector> missing;

    for (int i = 0; i < size_; i++)
        if (providers_[i] == -1)
            missing.append(points_[i]);

    return missing;
}

template class pointInterpolator<colocated>;
template class pointInterpolator<staggered>;

}

}

}
