#ifndef midPointReconstructionScheme_H
#define midPointReconstructionScheme_H

#include "reconstructionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class midPointReconstructionScheme
:
    public reconstructionScheme<Type>
{
public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from dictionary and mesh
    midPointReconstructionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointReconstructionScheme(const fvMesh&);

    // Reconstruction of a staggered quantity into a colocated one
    virtual tmp
    <
        meshField
        <
            typename reconstructionScheme<Type>::ReconType,
            colocated
        >
    >
    reconstruct
    (
        const meshField<Type,staggered>&
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointReconstructionScheme.C"
#endif

#endif
