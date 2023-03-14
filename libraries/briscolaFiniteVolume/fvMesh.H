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
    public regIOobject
{
    // Pointer to the mesh. The mesh cannot be a base class, because it is
    // runtime selectable.
    autoPtr<mesh> mshPtr_;

    // Scheme dictionary
    IOdictionary schemeDict_;

    // Solver dictionary
    IOdictionary solverDict_;

    // Colocated mesh metrics
    autoPtr<fvMeshMetrics<colocated>> colocatedMetrics_;

    // Staggered mesh metrics
    autoPtr<fvMeshMetrics<staggered>> staggeredMetrics_;

    // Disallow default bitwise assignment
    void operator=(const fvMesh&);

public:

    // Declare name of the class and its debug switch
    ClassName("fvMesh");

    // Construct from IOdictionary
    fvMesh(const IOdictionary&, const Time&);

    // Copy constructor
    fvMesh(const fvMesh&);

    // Destructor
    ~fvMesh();

    // Return access to the mesh
    inline mesh& msh()
    {
        return mshPtr_();
    }

    // Return const access to the mesh
    inline const mesh& msh() const
    {
        return mshPtr_();
    }

    // Const access to face centers
    template<class MeshType>
    const fvMeshMetrics<MeshType>& metrics() const;

    // Scheme dictionary
    inline const dictionary& schemeDict() const { return schemeDict_; }

    // Solver dictionary
    inline const dictionary& solverDict() const { return solverDict_; }

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Member functions that pass directly to the mesh object

        // Access to the part level
        inline partLevel& operator[](const label i)
        {
            return mshPtr_->operator[](i);
        }

        // Const access to the part level
        inline const partLevel& operator[](const label i) const
        {
            return mshPtr_->operator[](i);
        }

        // Number of levels
        inline label size() const
        {
            return mshPtr_->size();
        }

        // Const access to part patched
        inline const PtrList<partPatch>& partPatches() const
        {
            return mshPtr_->partPatches();
        }

        // Structured mesh bool
        inline bool structured() const
        {
            return mshPtr_->structured();
        }

        // Return rectilinear bools
        inline labelVector rectilinear() const
        {
            return mshPtr_->rectilinear();
        }

        // Return uniform bools
        inline labelVector uniform() const
        {
            return mshPtr_->uniform();
        }
};

}

}

}

#endif
