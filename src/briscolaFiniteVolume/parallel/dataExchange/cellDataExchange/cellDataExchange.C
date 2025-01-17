#include "cellDataExchange.H"
#include "domainBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(cellDataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(cellDataExchange<staggered>, 0);

template<class MeshType>
cellDataExchange<MeshType>::cellDataExchange
(
    const List<labelVector>& indices,
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    dataExchange<MeshType>(fvMsh,l,d)
{
    init(indices);
}

template<class MeshType>
cellDataExchange<MeshType>::cellDataExchange
(
    const cellDataExchange& e
)
:
    dataExchange<MeshType>(e),
    indices_(e.indices_),
    map_(e.map_),
    neighbors_(e.neighbors_),
    sendCells_(e.sendCells_),
    recvCells_(e.recvCells_)
{}

template<class MeshType>
cellDataExchange<MeshType>::~cellDataExchange()
{}

template<class MeshType>
void cellDataExchange<MeshType>::init(const List<labelVector>& indices)
{
    indices_ = indices;

    // Nothing to prepare when this is a serial run

    if (!Pstream::parRun())
        return;

    const faceLabel I = this->fvMsh_.template I<MeshType>(this->l_,this->d_);

    const labelVector U = I.upper();
    const labelVector L = I.lower();

    const mesh& msh = this->fvMsh_.msh();

    // Set neighboring processor numbers, offsets and transformations. The first
    // one is us.

    neighbors_.clear();
    neighbors_.setSize(msh.boundaries().size()+1,-1);
    neighbors_[0] = Pstream::myProcNo();

    List<labelVector> neighborOffsets(neighbors_.size(), zeroXYZ);
    List<labelTensor> neighborTs(neighbors_.size(), eye);

    forAll(msh.boundaries(), i)
    {
        if (msh.boundaries()[i].castable<parallelBoundary>())
        {
            neighbors_[i+1] = msh.boundaries()[i].neighborProcNum();
            neighborOffsets[i+1] = msh.boundaries()[i].offset();
            neighborTs[i+1] = msh.boundaries()[i].T();
        }
    }

    // Collect the required cells and store per neighbor processor

    recvCells_.clear();
    recvCells_.setSize(msh.boundaries().size()+1);

    map_.clear();
    map_.setSize(indices.size());

    forAll(indices, i)
    {
        labelVector index = indices[i];

        const labelVector closest = labelVector
        (
            Foam::max(Foam::min(index.x(), U.x()), L.x()-1),
            Foam::max(Foam::min(index.y(), U.y()), L.y()-1),
            Foam::max(Foam::min(index.z(), U.z()), L.z()-1)
        );

        const labelVector offset = index - closest;

        if (offset == zeroXYZ)
        {
            // Local point

            recvCells_[0].append(index);
            map_[i] = Tuple2<label,label>(0,recvCells_[0].size()-1);
        }
        else
        {
            const labelVector bo =
                briscola::cmptMin(briscola::cmptMax(offset, -unitXYZ), unitXYZ);

            const label degree = cmptSum(cmptMag(bo));

            if (!this->fvMsh_.structured() && degree > 1)
                FatalErrorInFunction
                    << "Can only obtain cell data across faces "
                    << "on unstructured meshes" << endl
                    << abort(FatalError);

            // Get the boundary to the other processor

            label bNum;
            for (bNum = 1; bNum < neighborOffsets.size(); bNum++)
                if (neighborOffsets[bNum] == bo)
                    break;

            // Error when not found

            if (bNum == neighborOffsets.size())
                FatalErrorInFunction
                    << "Could not find boundary corresponding to offset "
                    << bo << endl << abort(FatalError);

            // Error when we're outside the domain

            if
            (
                msh.boundaries()[bNum-1].castable<domainBoundary>()
             || msh.boundaries()[bNum-1].castable<emptyBoundary>()
            )
                FatalErrorInFunction
                    << "Cannot exchange cell data across a domain boundary "
                    << "at offset " << bo << endl << abort(FatalError);

            // Store trimmed index

            for (int j = 0; j < 3; j++)
                if (bo[j] != 0)
                    index[j] -= bo[j] > 0 ? (U[j]-1) :L[j];

            recvCells_[bNum].append(index);

            map_[i] =
                Tuple2<label,label>
                (
                    bNum,
                    recvCells_[bNum].size()-1
                );
        }
    }

    // Exchange processor numbers, sizes and boundary offsets

    List<List<FixedList<label,5>>> globalData(Pstream::nProcs());

    forAll(recvCells_, i)
    {
        if (recvCells_[i].size() > 0)
        {
            FixedList<label,5> data(0);

            data[0] = neighbors_[i];
            data[1] = recvCells_[i].size();

            // Compute the boundary offset that the neighbor should have

            if (i > 0)
            {
                const labelVector bo = - (neighborTs[i].T() & neighborOffsets[i]);

                for (int j = 0; j < 3; j++)
                    data[j+2] = bo[j];
            }

            globalData[Pstream::myProcNo()].append(data);
        }
    }

    Pstream::gatherList(globalData);
    Pstream::scatterList(globalData);

    // Collect

    sendCells_.clear();
    sendCells_.setSize(recvCells_.size());

    forAll(globalData, proc)
    forAll(globalData[proc], i)
    {
        if (globalData[proc][i][0] == Pstream::myProcNo())
        {
            // Requested boundary offset

            const labelVector bo
            (
                globalData[proc][i][2],
                globalData[proc][i][3],
                globalData[proc][i][4]
            );

            // Find the matching neighbor

            bool found = false;

            forAll(neighbors_, j)
            if (neighbors_[j] == proc && bo == neighborOffsets[j])
            {
                sendCells_[j].setSize
                (
                    globalData[proc][i][1]
                );

                found = true;
                break;
            }

            if (!found)
                FatalErrorInFunction
                    << "Could not find matching neighbor for boundary offset "
                    << bo << endl << abort(FatalError);
        }
    }

    forAll(sendCells_, i)
        if (sendCells_[i].size() > 0)
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                neighbors_[i],
                reinterpret_cast<char*>(sendCells_[i].begin()),
                sendCells_[i].byteSize()
            );

    forAll(recvCells_, i)
        if (recvCells_[i].size() > 0)
            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                neighbors_[i],
                reinterpret_cast<char*>(recvCells_[i].begin()),
                recvCells_[i].byteSize()
            );

    UOPstream::waitRequests();

    // Displace and transform send cell indices. Not needed for local cells.

    forAll(sendCells_, i)
    if (i != 0 && sendCells_[i].size() > 0)
    {
        const labelTensor T = neighborTs[i];
        const faceLabel J = pTransform<faceLabel>(T,I);

        const labelVector L = briscola::cmptMin(J.lower(), J.upper());
        const labelVector U = briscola::cmptMax(J.lower(), J.upper());

        const labelVector bo = (T & neighborOffsets[i]);

        forAll(sendCells_[i], j)
        {
            labelVector& index = sendCells_[i][j];

            for (int k = 0; k < 3; k++)
                if (bo[k] != 0)
                    index[k] += index[k] > 0 ? (L[k]-1) : U[k];

            // Transform back

            index = (T.T() & index);

            for (int k = 0; k < 3; k++)
                index[k] += index[k] < 0 ? U[k] : 0;
        }
    }
}

template<class MeshType>
template<class Type>
List<Type> cellDataExchange<MeshType>::dataFunc
(
    const meshField<Type,MeshType>& field
) const
{
    List<Type> D(indices_.size());

    if (!Pstream::parRun())
    {
        forAll(D, i)
            D[i] = field(this->l_, this->d_, indices_[i]);

        return D;
    }
    else
    {
        List<List<Type>> data(recvCells_.size());

        forAll(recvCells_, i)
        if (recvCells_[i].size() > 0)
        {
            data[i].setSize(recvCells_[i].size());

            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                neighbors_[i],
                reinterpret_cast<char*>(data[i].begin()),
                data[i].byteSize()
            );
        }

        forAll(sendCells_, i)
        if (sendCells_[i].size() > 0)
        {
            List<Type> data(sendCells_[i].size());

            forAll(data, j)
                data[j] = field(this->l_, this->d_, sendCells_[i][j]);

            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                neighbors_[i],
                reinterpret_cast<char*>(data.begin()),
                data.byteSize()
            );
        }

        UOPstream::waitRequests();

        forAll(D, i)
            D[i] = data[map_[i].first()][map_[i].second()];

        return D;
    }
}

// Instantiate class and data functions

template class cellDataExchange<colocated>;
template class cellDataExchange<staggered>;

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
template List<TYPE> cellDataExchange<MESHTYPE>::dataFunc                       \
(                                                                              \
    const meshField<TYPE,MESHTYPE>&                                            \
) const;

INSTANTIATE(scalar,colocated)
INSTANTIATE(vector,colocated)
INSTANTIATE(tensor,colocated)
INSTANTIATE(sphericalTensor,colocated)
INSTANTIATE(symmTensor,colocated)
INSTANTIATE(diagTensor,colocated)
INSTANTIATE(faceScalar,colocated)
INSTANTIATE(faceVector,colocated)

INSTANTIATE(scalar,staggered)
INSTANTIATE(vector,staggered)
INSTANTIATE(tensor,staggered)
INSTANTIATE(sphericalTensor,staggered)
INSTANTIATE(symmTensor,staggered)
INSTANTIATE(diagTensor,staggered)
INSTANTIATE(faceScalar,staggered)
INSTANTIATE(faceVector,staggered)

#undef INSTANTIATE

}

}

}
