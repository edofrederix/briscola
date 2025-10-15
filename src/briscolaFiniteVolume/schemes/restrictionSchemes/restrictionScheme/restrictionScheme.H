#ifndef restrictionScheme_H
#define restrictionScheme_H

#include "scheme.H"
#include "faceFields.H"
#include "meshFields.H"
#include "meshLevels.H"
#include "meshDirections.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class restrictionScheme
:
    public scheme
{
protected:

    // Raise error for restriction schemes that are incompatible with scaling
    void errorNoScaling(const bool scale)
    {
        if (scale)
            FatalErrorInFunction
                << "Restriction scheme type " << this->typeName
                << " cannot be used as a scaled restriction scheme"
                << endl << abort(FatalError);
    }

public:

    // Runtime type information
    TypeName("restrictionScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        restrictionScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    restrictionScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    restrictionScheme(const restrictionScheme<Type,MeshType>&);

    // Selector
    static autoPtr<restrictionScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~restrictionScheme();

    // Default restriction scheme for this mesh and type
    static const word defaultScheme;

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level, with option to scale by cell volume
    virtual void restrict
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const bool = false
    ) = 0;

    // Restrict from tmp and clear tmp
    virtual void restrict
    (
        meshDirection<Type,MeshType>& res,
        const tmp<meshDirection<Type,MeshType>>& tf,
        const bool scale = false
    );

    // Restrict level
    virtual void restrict
    (
        meshLevel<Type,MeshType>& res,
        const meshLevel<Type,MeshType>& f,
        const bool scale = false
    );

    // Restrict tmp level and clear tmp
    virtual void restrict
    (
        meshLevel<Type,MeshType>& res,
        const tmp<meshLevel<Type,MeshType>>& tf,
        const bool scale = false
    );

    // Restrict field from fine to coarse
    virtual void restrict
    (
        meshField<Type,MeshType>& res,
        const meshField<Type,MeshType>& f,
        const bool scale = false
    );

    // Restrict tmp field from fine to coarse and clear tmp
    virtual void restrict
    (
        meshField<Type,MeshType>& res,
        const tmp<meshField<Type,MeshType>>& tf,
        const bool scale = false
    );

    // Restrict field from itself from fine to coarse
    virtual void restrict
    (
        meshField<Type,MeshType>& res,
        const bool scale = false
    );
};

}

}

}

#ifdef NoRepository
    #include "restrictionScheme.C"
#endif

#endif
