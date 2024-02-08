#include "globalCellNumbers.H"

#include "meshField.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class MeshType>
globalCellNumbers<SType,MeshType>::globalCellNumbers
(
    const fvMesh& fvMsh,
    const label l
)
:
    fvMsh_(fvMsh),
    l_(l)
{
    data_.setSize(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        data_.set
        (
            d,
            new List<List<FixedList<label,SType::nComponents>>>
            (
                Pstream::nProcs()
            )
        );

        List<FixedList<label,SType::nComponents>>& data =
            data_[d][Pstream::myProcNo()];

        const meshDirection<label,MeshType>& numbers =
            fvMsh_.template metrics<MeshType>().globalCellNumbers()[l][d];

        data.setSize(cmptProduct(numbers.N()));

        int c = 0;
        forAllCells(numbers, i, j, k)
        {
            for(int s = 0; s < SType::nComponents; s++)
                data[c][s] =
                    numbers
                    (
                        labelVector(i,j,k)
                      + SType::componentOffsets[s]
                    );

            c++;
        }

        if (Pstream::parRun())
            Pstream::gatherList(data_[d]);
    }
}

template<class SType, class MeshType>
globalCellNumbers<SType,MeshType>::~globalCellNumbers()
{}

// Instantiate

template class globalCellNumbers<stencil,colocated>;
template class globalCellNumbers<symmStencil,colocated>;

template class globalCellNumbers<stencil,staggered>;
template class globalCellNumbers<symmStencil,staggered>;

}

}

}
