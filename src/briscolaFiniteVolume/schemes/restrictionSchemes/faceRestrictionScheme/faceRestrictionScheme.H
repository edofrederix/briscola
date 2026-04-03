#ifndef faceRestrictionScheme_H
#define faceRestrictionScheme_H

#include "scheme.H"
#include "faceFieldFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class faceRestrictionScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("faceRestrictionScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        faceRestrictionScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    faceRestrictionScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    faceRestrictionScheme(const faceRestrictionScheme<Type,MeshType>&);

    // Selector from type name
    static autoPtr<faceRestrictionScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~faceRestrictionScheme();

    // Default restriction scheme for this mesh and type
    static const word defaultScheme;

    // In-place restriction from a fine face level l to an adjacent coarse face
    // level l+1, for given direction d
    virtual void restrict
    (
        faceField<Type,MeshType>&,
        const label l
    ) = 0;

    // Restrict field from itself from fine to coarse
    virtual void restrict
    (
        faceField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "faceRestrictionScheme.C"
#endif

#endif
