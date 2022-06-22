#ifndef scheme_H
#define scheme_H

#include "fvMesh.H"
#include "dictionary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class scheme
{
    // Copy of the scheme dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

public:

    // Runtime type information
    TypeName("scheme");

    // Construct from dictionary
    scheme(const dictionary&, const fvMesh& fvMsh);

    // Copy constructor
    scheme(const scheme&);

    // Destructor
    virtual ~scheme();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }
};

}

}

}

#endif
