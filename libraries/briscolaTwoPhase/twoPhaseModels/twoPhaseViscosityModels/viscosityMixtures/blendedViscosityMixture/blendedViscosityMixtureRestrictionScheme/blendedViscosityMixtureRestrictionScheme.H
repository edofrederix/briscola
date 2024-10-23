#ifndef blendedViscosityMixtureRestrictionScheme_H
#define blendedViscosityMixtureRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class blendedViscosityMixtureRestrictionScheme
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
        const scalar xHat(Foam::min(Foam::max(x, 1e-12), 1-1e-12));

        return
            (mu2_ - mu1_)/2.0
          * (Foam::tanh(C_*Foam::atanh(2.0*xHat - 1.0)) - 1.0)
          + mu2_;
    }

    // The inverse blending function
    inline scalar inv(const scalar y) const
    {
        const scalar xHat
        (
            Foam::min(Foam::max((y - mu2_)/(mu1_ - mu2_), 1e-12), 1-1e-12)
        );

        return
            0.5
          * (
                Foam::tanh
                (
                    Foam::atanh
                    (
                        1.0 - 2.0*xHat
                    )/C_
                )
              + 1.0
            );
    }

public:

    // Runtime type information
    TypeName("blendedViscosityMixture");

    // Construct from stream and mesh
    blendedViscosityMixtureRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosityMixtureRestrictionScheme.C"
#endif

#endif
