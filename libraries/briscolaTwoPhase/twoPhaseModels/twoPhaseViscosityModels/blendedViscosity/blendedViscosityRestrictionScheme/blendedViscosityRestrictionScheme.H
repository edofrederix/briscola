#ifndef blendedViscosityRestrictionScheme_H
#define blendedViscosityRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class blendedViscosityRestrictionScheme
:
    public restrictionScheme<lowerFaceScalar,MeshType>
{
protected:

    // Blending function coefficients taken from the viscosity model
    scalar C_;
    scalar mu1_;
    scalar mu2_;

    // The blending function
    inline scalar blend(const scalar x) const
    {
        return
            (mu2_ - mu1_)/2.0
          * (Foam::tanh(C_*Foam::atanh(2.0*x - 1.0)) - 1.0)
          + mu2_;
    }

    // The inverse blending function
    inline scalar inv(const scalar y) const
    {
        return
            0.5
          * (
                Foam::tanh
                (
                    Foam::atanh
                    (
                        1.0 - 2.0*(y - mu2_)/(mu1_ - mu2_)
                    )/C_
                )
              + 1.0
            );
    }

public:

    // Runtime type information
    TypeName("blendedViscosity");

    // Construct from stream and mesh
    blendedViscosityRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosityRestrictionScheme.C"
#endif

#endif
