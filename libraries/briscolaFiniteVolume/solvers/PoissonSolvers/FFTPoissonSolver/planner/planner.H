#ifndef planner_H
#define planner_H

#include "List.H"
#include "scalarList.H"
#include "vector.H"
#include "block.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"
#include "tridiagonalSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace FFT
{

class planner
{
    // Initial decomposition
    labelVector I_;

    // Mesh dimensions
    labelVector N_;

    // Initial decomposition type:
    // 0: Not decomposed
    // 1: 3D
    // 2: x-pencil
    // 3: y-pencil
    // 4: z-pencil
    // 5: xy-slab
    // 6: xz-slab
    // 7: yz-slab
    label decompType_;

    // Solve direction
    label solveDir_;

public:

    // Construct
    planner
    (
        const fvMesh& fvMsh
    );

    // Destructor
    virtual ~planner();

    inline label decompType() const { return decompType_; }
    inline label solveDir() const { return solveDir_; }

};

}

}

}

}

#endif
