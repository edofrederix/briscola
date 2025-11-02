#include "createBriscolaMesh.H"

if (!fvMsh.msh().castable<rectilinearMesh>())
    FatalErrorInFunction
        << "Cannot create staggered mesh. Mesh must be rectilinear."
        << abort(FatalError);
