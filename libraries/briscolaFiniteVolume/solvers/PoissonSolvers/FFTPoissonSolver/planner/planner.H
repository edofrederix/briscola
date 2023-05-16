#ifndef planner_H
#define planner_H

#include "List.H"
#include "scalarList.H"
#include "vector.H"
#include "block.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"
#include "FourierTransforms.H"
#include "tridiagonalSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
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

    // Mesh uniformity per direction (0 = false, 1 = true)
    labelVector meshUniform_;

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

    inline label decompType() { return decompType_; }
    inline labelVector meshUniform() { return meshUniform_; }
    inline label solveDir() { return solveDir_; }

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
