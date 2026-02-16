#ifndef mesh_H
#define mesh_H

#include "decomposition.H"
#include "geometry.H"
#include "level.H"
#include "faceLabel.H"
#include "edgeLabel.H"
#include "vertexLabel.H"

namespace Foam
{

namespace briscola
{

class mesh
:
    public geometry,
    public FastPtrList<level>
{
    // Bool to check if mesh is structured or unstructured. If the brick
    // topology is structured, so is the mesh.
    bool structured_;

    // Check if mesh directions are rectilinear
    labelVector rectilinear_;

    // Check if mesh directions are uniform
    labelVector uniform_;

    // Bounding box
    faceScalar boundingBox_;

    // Generate mesh levels
    void generateLevels();

public:

    // Runtime type information
    TypeName("mesh");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        mesh,
        dictionary,
        (autoPtr<mesh>& mshPtr),
        (mshPtr)
    );

    // Construct from IOdictionary
    explicit mesh(const IOdictionary&);

    // Copy constructor
    mesh(const mesh& msh);

    // Transfer constructor and reset arg
    mesh(autoPtr<mesh>&);

    // Clone mesh
    virtual autoPtr<mesh> clone() const
    {
        return autoPtr<mesh>(new mesh(*this));
    }

    // Selector
    static autoPtr<mesh> New(const IOdictionary&);

    // Destructor
    virtual ~mesh();

    // Check if cast to derived type is possible
    template<class MeshType>
    inline bool castable() const
    {
        return dynamic_cast<const MeshType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class MeshType>
    inline const MeshType& cast() const
    {
        const MeshType* mshPtr = dynamic_cast<const MeshType*>(this);

        if (!mshPtr)
        {
            FatalErrorInFunction
                << "Mesh cannot be cast to " << MeshType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *mshPtr;
    }

    // Structured mesh bool
    inline bool structured() const { return structured_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // Bounding box
    inline faceScalar boundingBox() const { return boundingBox_; }

    // Decomposition of the first level
    inline const decomposition& decomp() const
    {
        return this->operator[](0).decomp();
    }

    // Boundaries of the first level
    inline const levelBoundaries& boundaries() const
    {
        return this->operator[](0).boundaries();
    }

    // Communicators of the first level
    inline const levelComms& comms() const
    {
        return this->operator[](0).comms();
    }

    // Communicator number of the first level
    inline label comm() const
    {
        return this->operator[](0).comms();
    }

    // Points of the first level
    inline const levelPoints& points() const
    {
        return this->operator[](0).points();
    }

    // Find the cell index of the given point at given mesh level. Return
    // -unitXYZ if it is not found. Tree search, but should be overloaded by
    // mesh specializations for better performance.
    virtual labelVector findCell(const vector& p, const label l = 0) const;

    // Find the cell indices of the given list of points at given mesh level
    virtual List<labelVector> findCells(const vectorList&, const label) const;
};

}

}

#endif
