#ifndef faceDotGradientScheme_H
#define faceDotGradientScheme_H

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
class faceDotGradientScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("faceDotGradientScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        faceDotGradientScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    faceDotGradientScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    faceDotGradientScheme(const faceDotGradientScheme<Type,MeshType>&);

    // Selector
    static autoPtr<faceDotGradientScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~faceDotGradientScheme();

    // Face-dot gradient multiplied by face area
    virtual tmp<faceField<Type,MeshType>> faceDotGrad
    (
        const meshField<Type,MeshType>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "faceDotGradientScheme.C"
#endif

#endif
