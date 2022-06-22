#ifndef fvMesh_H
#define fvMesh_H

#include "mesh.H"

#include "regIOobject.H"
#include "autoPtr.H"
#include "Time.H"

#include "fvMeshMetrics.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh
:
    public mesh,
    public regIOobject
{
    // Scheme dictionary
    IOdictionary schemeDict_;

    // Solver dictionary
    IOdictionary solverDict_;

    // Colocated mesh metrics
    autoPtr<fvMeshMetrics<colocated>> colocatedMetrics_;

    // Staggered mesh metrics
    autoPtr<fvMeshMetrics<staggered>> staggeredMetrics_;

public:

    // Declare name of the class and its debug switch
    ClassName("fvMesh");

    // Construct from IOdictionary
    explicit fvMesh(const IOdictionary&, const Time&);

    // Disallow default bitwise copy construction
    fvMesh(const fvMesh&);

    // Destructor
    ~fvMesh();

    // Disallow default bitwise assignment
    void operator=(const fvMesh&);

    // Const access to face centers
    template<class MeshType>
    const fvMeshMetrics<MeshType>& metrics() const;

    // Scheme dictionary
    inline const dictionary& schemeDict() const { return schemeDict_; }

    // Solver dictionary
    inline const dictionary& solverDict() const { return solverDict_; }

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }
};

}

}

}

#endif
