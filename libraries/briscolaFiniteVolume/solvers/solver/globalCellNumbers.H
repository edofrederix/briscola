#ifndef globalCellNumbers_H
#define globalCellNumbers_H

#include "PtrList.H"
#include "FixedList.H"
#include "List.H"
#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class MeshType>
class globalCellNumbers
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    const label l_;

    // Data for each direction, processor, cell and stencil component
    PtrList<List<List<FixedList<label,SType::nComponents>>>> data_;

public:

    // Construct from mesh and level number
    globalCellNumbers(const fvMesh&, const label);

    // Destructor
    virtual ~globalCellNumbers();

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Level number
    inline label levelNum() const { return l_; }

    // Return data
    inline const auto& data() const { return data_; }

    // Return directional data
    inline const auto& data(const label d) const { return data_[d]; }

    // Compute the total number of cells for given direction
    inline label numberOfCells(const label d) const
    {
        label N = 0;
        forAll(data_[d], proc)
            N += data_[d][proc].size();

        return N;
    }
};

}

}

}

#endif
