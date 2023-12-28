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

    const faceLabel I = this->fvMsh_.template I<MeshType>(this->l_,this->d_);

    const labelVector U = I.upper();
    const labelVector L = I.lower();

    const mesh& msh = this->fvMsh_.msh();

    // Set neighboring processor numbers and boundaries. The first one is us.

    neighbors_.clear();
    neighbors_.setSize(msh.boundaries().size()+1,-1);
    neighbors_[0] = Pstream::myProcNo();

    forAll(msh.boundaries(), i)
        if (msh.boundaries()[i].castable<parallelBoundary>())
            neighbors_[i+1] = msh.boundaries()[i].neighborProcNum();

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
            Foam::max(Foam::min(index.x(), U.x()-1), L.x()),
            Foam::max(Foam::min(index.y(), U.y()-1), L.y()),
            Foam::max(Foam::min(index.z(), U.z()-1), L.z())
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

            label bNum = -1;

            forAll(msh.boundaries(), bi)
            if (msh.boundaries()[bi].offset() == bo)
            {
                bNum = bi;
                break;
            }

            if (bNum == -1)
                FatalErrorInFunction
                    << "Could not find boundary corresponding to boundary "
                    << "offset " << bo << endl << abort(FatalError);

            const boundary& b = msh.boundaries()[bNum];

            if (b.castable<domainBoundary>() || b.castable<emptyBoundary>())
                FatalErrorInFunction
                    << "Cannot exchange cell data across a domain boundary"
                    << endl << abort(FatalError);

            // Store trimmed index

            for (int j = 0; j < 3; j++)
                if (bo[j] != 0)
                    index[j] = bo[j] > 0 ? index[j]-(U[j]-1) : index[j]-L[j];

            recvCells_[bNum+1].append(index);

            map_[i] =
                Tuple2<label,label>
                (
                    bNum+1,
                    recvCells_[bNum+1].size()-1
                );
        }
    }

    // Exchange processor numbers to talk to and data sizes

    List<List<Tuple2<label,label>>> globalData(Pstream::nProcs());

    forAll(recvCells_, i)
        if (recvCells_[i].size() > 0)
            globalData[Pstream::myProcNo()].append
            (
                Tuple2<label,label>(neighbors_[i], recvCells_[i].size())
            );

    Pstream::gatherList(globalData);
    Pstream::scatterList(globalData);

    // Collect

    sendCells_.clear();
    sendCells_.setSize(recvCells_.size());

    forAll(sendCells_, i)
        if (neighbors_[i] > -1)
            forAll(globalData[neighbors_[i]], j)
                if (globalData[neighbors_[i]][j].first() == Pstream::myProcNo())
                    sendCells_[i].setSize
                    (
                        globalData[neighbors_[i]][j].second()
                    );

    // Send indices to neighbors and receive from neighbors

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

    UPstream::waitRequests();

    // Displace and transform send cell indices. Not needed for local cells.

    forAll(sendCells_, i)
    if (i != 0 && sendCells_[i].size() > 0)
    {
        const boundary& b = msh.boundaries()[i-1];
        const labelTensor& T = b.T();

        faceLabel J = pTransform<faceLabel>(T,I);

        const labelVector L = briscola::cmptMin(J.lower(), J.upper());
        const labelVector U = briscola::cmptMax(J.lower(), J.upper());

        const labelVector bo = (T & b.offset());

        forAll(sendCells_[i], j)
        {
            labelVector& index = sendCells_[i][j];

            for (int k = 0; k < 3; k++)
                if (bo[k] != 0)
                    index[k] = index[k] > 0 ? (L[k]-1)+index[k] : U[k]+index[k];

            // Transform back

            index = (T.T() & index);

            for (int k = 0; k < 3; k++)
                index[k] = index[k] < 0 ? U[k]+index[k] : index[k];
        }
    }
}

template<class MeshType>
template<class Type>
List<Type> cellDataExchange<MeshType>::dataFunc
(
    const meshField<Type,MeshType>& field
)
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
            data[j] = field(this->l_,this->d_,sendCells_[i][j]);

        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            neighbors_[i],
            reinterpret_cast<char*>(data.begin()),
            data.byteSize()
        );
    }

    UPstream::waitRequests();

    List<Type> D(indices_.size());

    forAll(D, i)
        D[i] = data[map_[i].first()][map_[i].second()];

    return D;
}

// Instantiate class and data functions

template class cellDataExchange<colocated>;
template class cellDataExchange<staggered>;

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
template List<TYPE> cellDataExchange<MESHTYPE>::dataFunc                       \
(                                                                              \
    const meshField<TYPE,MESHTYPE>&                                            \
);

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
