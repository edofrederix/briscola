#ifndef coloFaceAreaWeightedRestrictionScheme_H
#define coloFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("faceAreaWeighted");

    // Construct from dictionary and mesh
    coloFaceAreaWeightedRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    coloFaceAreaWeightedRestrictionScheme(const fvMesh&);

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
    #include "coloFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
