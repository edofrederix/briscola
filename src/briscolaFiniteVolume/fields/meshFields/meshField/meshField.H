#ifndef meshField_H
#define meshField_H

#include "block.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "FastPtrList.H"

#include "linearSystemFwd.H"
#include "meshLevel.H"
#include "immersedBoundaryConditionFwd.H"
#include "cachedRefCount.H"
#include "direction.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

template<class Type, class MeshType>
class meshField
:
    public FastPtrList<meshLevel<Type,MeshType>>,
    public IOdictionary,
    public cachedRefCount
{
public:

    // Number of ghost cell layers
    static const label ghosts = meshDirection<Type,MeshType>::ghosts;

private:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Old time field pointer, explicitly set by setOldTime and never changed by
    // any other operation. It is also not copied/created during assignment or
    // construction from another field that has the old time field pointer set.
    meshField<Type,MeshType>* oldTimePtr_;

    // Immersed boundary condition
    FastPtrList<immersedBoundaryCondition<Type,MeshType>>
        immersedBoundaryConditions_;

    // Restriction scheme name
    word restrictionScheme_;

protected:

    // Allocate data and specify depth
    void allocate(const bool);

    // Set field pointers in levels, needed when data is transferred
    void setFieldPointers();

    // Transfer data from other mesh field
    void transfer(meshField<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("meshField");

    // Data type
    typedef Type dataType;

    // Mesh type
    typedef MeshType meshType;

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // List type
    typedef FastPtrList<meshLevel<Type,MeshType>> listType;

    // Main constructor: construct from name and mesh with optional read and
    // write options, object registration and depth

    meshField
    (
        const word& name,
        const fvMesh& fvMsh,
        const IOobject::readOption r = NO_READ,
        const IOobject::writeOption w = NO_WRITE,
        const bool registerObject = false,
        const bool deep = false
    );

    meshField
    (
        const IOobject& io,
        const fvMesh& fvMsh,
        const bool deep = false
    );

    // Copy constructors: construct from field or tmp of field with option to
    // clone or transfer the boundary conditions

    meshField
    (
        const meshField<Type,MeshType>& field,
        const bool registerObject = false
    );

    meshField
    (
        const word& name,
        const meshField<Type,MeshType>& field,
        const bool registerObject = false
    );

    meshField
    (
        const tmp<meshField<Type,MeshType>>& field,
        const bool registerObject = false
    );

    meshField
    (
        const word& name,
        const tmp<meshField<Type,MeshType>>& field,
        const bool registerObject = false
    );

    // Temporary field constructor
    static tmp<meshField<Type,MeshType>> New
    (
        const word& name,
        const fvMesh& fvMsh,
        const bool cached = false
    )
    {
        return
            fvMsh.template
            NewObject<meshField<Type,MeshType>>(name, cached);
    }

    // Clone
    inline autoPtr<meshField<Type,MeshType>> clone() const
    {
        return autoPtr<meshField<Type,MeshType>>
        (
            new meshField<Type,MeshType>(*this)
        );
    }

    // Destructor
    ~meshField();

    // Access to the levels
    inline meshLevel<Type,MeshType>& operator[](label l)
    {
        return listType::operator[](l);
    }

    // Const access to the levels
    inline const meshLevel<Type,MeshType>& operator[](label l) const
    {
        return listType::operator[](l);
    }

    // Access to the first level
    inline meshLevel<Type,MeshType>& level()
    {
        return listType::operator[](0);
    }

    // Const access to the first level
    inline const meshLevel<Type,MeshType>& level() const
    {
        return listType::operator[](0);
    }

    // Alias access function
    inline meshLevel<Type,MeshType>& level(const int l)
    {
        return listType::operator[](l);
    }

    // Alias access function
    inline const meshLevel<Type,MeshType>& level(const int l) const
    {
        return listType::operator[](l);
    }

    // Access to the first direction of the first level
    inline meshDirection<Type,MeshType>& direction()
    {
        return listType::operator[](0)[0];
    }

    // Const access to the first direction of the first level
    inline const meshDirection<Type,MeshType>& direction() const
    {
        return listType::operator[](0)[0];
    }

    // Access to the given direction of the first level
    inline meshDirection<Type,MeshType>& direction(const int d)
    {
        return listType::operator[](0)[d];
    }

    // Const access to the given direction of the first level
    inline const meshDirection<Type,MeshType>& direction(const int d) const
    {
        return listType::operator[](0)[d];
    }

    // Access to the given direction and level
    inline meshDirection<Type,MeshType>& direction
    (
        const int l,
        const int d
    )
    {
        return listType::operator[](l)[d];
    }

    // Const access to the given direction and level
    inline const meshDirection<Type,MeshType>& direction
    (
        const int l,
        const int d
    ) const
    {
        return listType::operator[](l)[d];
    }

    // Number of levels
    inline label size() const { return listType::size(); }

    // Alias to size
    inline label numberOfLevels() const { return size(); }

    // Number of directions
    inline label numberOfDirections() const
    {
        return MeshType::numberOfDirections;
    }

    // Return the finite volume mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the mesh
    inline const mesh& msh() const { return fvMsh_.msh(); }

    // Old time field pointer
    inline meshField<Type,MeshType>* oldTimePtr() const
    {
        if (oldTimePtr_ == nullptr)
            this->setOldTime();

        return oldTimePtr_;
    }

    // Access to the old time field
    inline meshField<Type,MeshType>& oldTime()
    {
        if (oldTimePtr_ == nullptr)
            this->setOldTime();

        return *oldTimePtr_;
    }

    inline const meshField<Type,MeshType>& oldTime() const
    {
        if (oldTimePtr_ == nullptr)
            const_cast<meshField<Type,MeshType>&>(*this).setOldTime();

        return *oldTimePtr_;
    }

    // Access to the boundary conditions
    inline FastPtrList<immersedBoundaryCondition<Type,MeshType>>&
    immersedBoundaryConditions()
    {
        return immersedBoundaryConditions_;
    }

    // Const access to the boundary conditions
    inline const FastPtrList<immersedBoundaryCondition<Type,MeshType>>&
    immersedBoundaryConditions() const
    {
        return immersedBoundaryConditions_;
    }

    // Access to the boundary conditions of the first level
    inline FastPtrList<boundaryCondition<Type,MeshType>>& boundaryConditions()
    {
        return level().boundaryConditions();
    }

    // Access to the boundary conditions of the first level
    inline const FastPtrList<boundaryCondition<Type,MeshType>>&
    boundaryConditions() const
    {
        return level().boundaryConditions();
    }

    // Add boundary conditions
    void addBoundaryConditions();

    // Add immersed boundary conditions
    void addImmersedBoundaryConditions();

    // Transfer boundary conditions from another field
    void transferBoundaryConditions(const meshField<Type,MeshType>&);

    // Transfer immersed boundary conditions from another field
    void transferImmersedBoundaryConditions(const meshField<Type,MeshType>&);

    // Correct boundary conditions
    void correctBoundaryConditions();

    // Prepare/evaluate/correct selected boundary conditions limited to a
    // boundary offset degree (defaults to 3)
    template<class Selector, int Degree = 3> void prepare();
    template<class Selector, int Degree = 3> void evaluate();
    template<class Selector, int Degree = 3> void correct();

    // Correct immersed boundary conditions
    void correctImmersedBoundaryConditions();

    // Correct agglomeration data
    void correctAggData();

    // Set the old time field. Works recursively, so if the old time field has
    // an old time field too, it will be set as well.
    void setOldTime();

    // Make the field deep
    void makeDeep();

    // Allow a const field to be made deep
    void makeDeep() const
    {
        const_cast<meshField<Type,MeshType>*>(this)->makeDeep();
    }

    // Make the field shallow
    void makeShallow();

    // Allow a const field to be made shallow
    void makeShallow() const
    {
        const_cast<meshField<Type,MeshType>*>(this)->makeShallow();
    }

    // Make the field deep or shallow
    inline void make(const bool deep)
    {
        if (deep)
            makeDeep();
        else
            makeShallow();
    }

    // Allow a const field to be made deep or shallow
    inline void make(const bool deep) const
    {
        const_cast<meshField<Type,MeshType>*>(this)->make(deep);
    }

    // The field is deep if it contains more than one level
    inline bool deep() const
    {
        return listType::size() > 1;
    }

    // Shallow if not deep
    inline bool shallow() const
    {
        return !deep();
    }

    // Set restriction scheme
    void setRestrictionScheme(const word scheme)
    {
        restrictionScheme_ = scheme;
    }

    // Return restriction scheme name
    const word& restrictionScheme() const { return restrictionScheme_; }

    // Write data (not supported)
    bool writeData(Ostream& os) const { return os; }

    // Return a component
    tmp<meshField<cmptType,MeshType>> component(const label) const;

    // Set a component from a directional list of values
    void replace(const label, const List<cmptType>&);

    // Set a component from a component field
    void replace(const label, const meshField<cmptType,MeshType>&);
    void replace(const label, const tmp<meshField<cmptType,MeshType>>&);

    // Set maximum/minimum of field and given value
    void max(const Type&);
    void min(const Type&);

    // Member operators

    void operator=(const meshField<Type,MeshType>&);
    void operator=(const tmp<meshField<Type,MeshType>>&);
    void operator=(const Type&);
    void operator=(const List<Type>&);
    void operator=(const zero);

    void operator+=(const meshField<Type,MeshType>&);
    void operator+=(const tmp<meshField<Type,MeshType>>&);

    void operator-=(const meshField<Type,MeshType>&);
    void operator-=(const tmp<meshField<Type,MeshType>>&);

    void operator*=(const meshField<scalar,MeshType>&);
    void operator*=(const tmp<meshField<scalar,MeshType>>&);

    void operator/=(const meshField<scalar,MeshType>&);
    void operator/=(const tmp<meshField<scalar,MeshType>>&);

    void operator+=(const Type&);
    void operator+=(const List<Type>&);

    void operator-=(const Type&);
    void operator-=(const List<Type>&);

    void operator*=(const scalar&);
    void operator*=(const scalarList&);

    void operator/=(const scalar&);
    void operator/=(const scalarList&);

    // Templated functions are needed for types which can be cast to another
    // type

    #define TYPE2FIELD meshField<Type2,MeshType>

    template<class Type2> void operator=(const TYPE2FIELD&);
    template<class Type2> void operator=(const tmp<TYPE2FIELD>&);

    template<class Type2> void operator+=(const TYPE2FIELD&);
    template<class Type2> void operator+=(const tmp<TYPE2FIELD>&);

    template<class Type2> void operator-=(const TYPE2FIELD&);
    template<class Type2> void operator-=(const tmp<TYPE2FIELD>&);

    template<class Type2> void operator*=(const TYPE2FIELD&);
    template<class Type2> void operator*=(const tmp<TYPE2FIELD>&);

    template<class Type2> void operator/=(const TYPE2FIELD&);
    template<class Type2> void operator/=(const tmp<TYPE2FIELD>&);

    #undef TYPE2FIELD

    template<class Type2> void operator=(const Type2&);
    template<class Type2> void operator=(const List<Type2>&);

    template<class Type2> void operator+=(const Type2&);
    template<class Type2> void operator+=(const List<Type2>&);

    template<class Type2> void operator-=(const Type2&);
    template<class Type2> void operator-=(const List<Type2>&);

    template<class Type2> void operator*=(const Type2&);
    template<class Type2> void operator*=(const List<Type2>&);

    template<class Type2> void operator/=(const Type2&);
    template<class Type2> void operator/=(const List<Type2>&);

    // Befriend linear system
    template<class SType, class, class>
    friend class linearSystem;

    // Befriend mesh level
    template<class, class>
    friend class meshLevel;

    // Befriend mesh field set
    template<class, class, ::Foam::direction>
    friend class meshFieldSet;

    #include "meshLevelFirstDirectionFunctions.H"
    #include "meshLevelDirectionFunctions.H"
    #include "meshFieldDirectionFunctions.H"
};

}

}

}

#include "meshFieldFunctions.H"
#include "meshFieldStencilFunctions.H"

#ifdef NoRepository
    #include "meshField.C"
#endif

#endif
