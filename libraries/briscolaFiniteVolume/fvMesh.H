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

    // Colocated internal cell start and end indices
    List<List<faceLabel>> Ic_;

    // Staggered internal cell start and end indices
    List<List<faceLabel>> Is_;

    // Colocated mesh metrics
    autoPtr<fvMeshMetrics<colocated>> colocatedMetrics_;

    // Staggered mesh metrics
    autoPtr<fvMeshMetrics<staggered>> staggeredMetrics_;

    // Set the internal cell start and end indices
    template<class MeshType>
    void setInternalCells();

    // Access to the internal cell start and end indices
    template<class MeshType>
    faceLabel& I(const label l, const label d);

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

    // Const access to the internal cell start and end indices

    template<class MeshType>
    const faceLabel& I(const label l, const label d) const;

    template<class MeshType>
    inline const faceLabel& I(const label d) const
    {
        return I<MeshType>(0,d);
    }

    template<class MeshType>
    inline const faceLabel& I() const
    {
        return I<MeshType>(0,0);
    }

    // Internal cell block size

    template<class MeshType>
    inline labelVector N(const label l, const label d) const
    {
        const faceLabel& I = this->I<MeshType>(l,d);
        return I.upper()-I.lower();
    }

    template<class MeshType>
    inline labelVector N(const label d) const
    {
        return N<MeshType>(0,d);
    }

    template<class MeshType>
    inline labelVector N() const
    {
        return N<MeshType>(0,0);
    }

    // Return the start point of the internal cell boundary at given offset

    template<class MeshType>
    inline labelVector S
    (
        const label l,
        const label d,
        const labelVector& bo
    ) const
    {
        const faceLabel& I = this->I<MeshType>(l,d);

        labelVector S(I.lower());

        for (label d = 0; d < 3; d++)
            if (bo[d] == 1)
                S[d] = I.upper()[d]-1;

        return S;
    }

    template<class MeshType>
    inline labelVector S(const label d, const labelVector& bo) const
    {
        return S<MeshType>(0,d,bo);
    }

    template<class MeshType>
    inline labelVector S(const labelVector& bo) const
    {
        return S<MeshType>(0,0,bo);
    }

    // Return the end point of the internal cell boundary at given offset

    template<class MeshType>
    inline labelVector E
    (
        const label l,
        const label d,
        const labelVector& bo
    ) const
    {
        const faceLabel& I = this->I<MeshType>(l,d);

        labelVector E(I.upper());

        for (label d = 0; d < 3; d++)
            if (bo[d] == -1)
                E[d] = I.lower()[d]+1;

        return E;
    }

    template<class MeshType>
    inline labelVector E(const label d, const labelVector& bo) const
    {
        return E<MeshType>(0,d,bo);
    }

    template<class MeshType>
    inline labelVector E(const labelVector& bo) const
    {
        return E<MeshType>(0,0,bo);
    }

    // Return the shape of the internal cell boundary at given offset

    template<class MeshType>
    inline labelVector N
    (
        const label l,
        const label d,
        const labelVector& bo
    ) const
    {
        return E<MeshType>(l,d,bo) - S<MeshType>(l,d,bo);
    }

    template<class MeshType>
    inline labelVector N(const label d, const labelVector& bo) const
    {
        return N<MeshType>(0,d,bo);
    }

    template<class MeshType>
    inline labelVector N(const labelVector& bo) const
    {
        return N<MeshType>(0,0,bo);
    }

    // Const access to face centers
    template<class MeshType>
    const fvMeshMetrics<MeshType>& metrics() const;

    // Cell find function templated for the mesh type. Uses the cell find
    // function of the base mesh and always acts on the first level.
    template<class MeshType>
    labelVector findCell
    (
        const vector& p,
        const label l = 0,
        const label d = 0
    ) const;

    // Helper to the find cell function
    template<class MeshType>
    List<labelVector> findCells
    (
        const vectorList&,
        const label l = 0,
        const label d = 0
    ) const;

    // Scheme dictionary
    inline const dictionary& schemeDict() const { return schemeDict_; }

    // Solver dictionary
    inline const dictionary& solverDict() const { return solverDict_; }

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Member functions that pass directly to the mesh object

        // Access to the mesh part
        inline part& operator[](const label i)
        {
            return mshPtr_->operator[](i);
        }

        // Const access to the mesh part
        inline const part& operator[](const label i) const
        {
            return mshPtr_->operator[](i);
        }

        // Number of levels
        inline label size() const
        {
            return mshPtr_->size();
        }

        // Const access to boundaries
        inline const PtrList<boundary>& boundaries() const
        {
            return mshPtr_->boundaries();
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
