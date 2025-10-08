#ifndef blendedViscosityMixtureRestrictionSchemeBase_H
#define blendedViscosityMixtureRestrictionSchemeBase_H

#include "restrictionScheme.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class blendedViscosityMixtureRestrictionSchemeBase
:
    public restrictionScheme<faceScalar,MeshType>
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
            Foam::min
            (
                Foam::max
                (
                    (y - mu2_)/stabilise(mu1_ - mu2_, 1e-12),
                    1e-12
                ),
                1-1e-12
            )
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

    // Construct from stream and mesh
    blendedViscosityMixtureRestrictionSchemeBase
    (
        const fvMesh& fvMsh,
        Istream& is = scheme::nullStream
    )
    :
        restrictionScheme<faceScalar,MeshType>(fvMsh, is)
    {
        const twoPhaseModel& model =
            fvMsh.db().lookupObject<twoPhaseModel>("briscolaTwoPhaseDict");

        C_ = model.dict().lookupOrDefault<scalar>("C", 8);

        // This currently assumes a Newtonian fluid
        mu1_ = readScalar(model.dict().lookup("mu1"));
        mu2_ = readScalar(model.dict().lookup("mu2"));
    }
};

}

}

}

#endif
