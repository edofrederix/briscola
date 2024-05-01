#ifndef coloHarmonicFaceAreaWeightedRestrictionScheme_H
#define coloHarmonicFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloHarmonicFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

    // Construct from dictionary and mesh
    coloHarmonicFaceAreaWeightedRestrictionScheme
    (
        const dictionary&,
        const fvMesh&
    );

    // Construct from mesh only, with empty dictionary
    coloHarmonicFaceAreaWeightedRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<Type,colocated>&,
        const meshDirection<Type,colocated>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "coloHarmonicFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
