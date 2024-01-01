#include "globalRowData.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
globalRowData<SType,Type,MeshType>::globalRowData
(
    const linearSystem<SType,Type,MeshType>& sys,
    const label l
)
:
    fvMsh_(sys.fvMsh()),
    l_(l)
{
    data_.setSize(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        data_.set
        (
            d,
            new List<List<Row<SType,Type>>>
            (
                Pstream::nProcs()
            )
        );

        List<Row<SType,Type>>& data = data_[d][Pstream::myProcNo()];

        data.setSize(cmptProduct(fvMsh_.template N<MeshType>(l,d)));

        const meshDirection<SType,MeshType>& A = sys.A()[l][d];
        const meshDirection<Type,MeshType>& b = sys.b()[l][d];

        int c = 0;
        forAllCells(A, i, j, k)
            data[c++] = Row<SType,Type>(A(i,j,k),b(i,j,k));

        Pstream::gatherList(data_[d]);
    }
}

template<class SType, class Type, class MeshType>
globalRowData<SType,Type,MeshType>::~globalRowData()
{}

// Instantiate

template class globalRowData<stencil,scalar,colocated>;
template class globalRowData<symmStencil,scalar,colocated>;
template class globalRowData<stencil,vector,colocated>;
template class globalRowData<symmStencil,vector,colocated>;

template class globalRowData<stencil,scalar,staggered>;
template class globalRowData<symmStencil,scalar,staggered>;
template class globalRowData<stencil,vector,staggered>;
template class globalRowData<symmStencil,vector,staggered>;

}

}

}
