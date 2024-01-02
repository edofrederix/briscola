#ifndef globalRowData_H
#define globalRowData_H

#include "PtrList.H"
#include "List.H"
#include "fvMesh.H"
#include "linearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class globalRowData
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    const label l_;

    // Row data for each direction, processor and cell
    PtrList<List<List<SType>>> data_;

public:

    // Construct from linear system and level number
    globalRowData(const linearSystem<SType,Type,MeshType>&, const label);

    // Destructor
    virtual ~globalRowData();

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Level number
    inline label levelNum() const { return l_; }

    // Return data
    inline const auto& data() const { return data_; }

    // Return directional data
    inline const auto& data(const label d) const { return data_[d]; }
};

}

}

}

#endif
