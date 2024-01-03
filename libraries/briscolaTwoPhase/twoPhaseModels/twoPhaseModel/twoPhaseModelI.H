#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

inline colocatedScalarField& twoPhaseModel::alpha()
{
    return alpha_;
}

inline const colocatedScalarField& twoPhaseModel::alpha() const
{
    return alpha_;
}

template<class MeshType>
inline tmp<meshField<scalar,MeshType>> twoPhaseModel::minRho() const
{
    tmp<meshField<scalar,MeshType>> tMinRho
    (
        min(this->rho1<MeshType>(), this->rho2<MeshType>())
    );
    tMinRho->setRestrictionScheme("volumeWeighted");
    return tMinRho;
}

template<class MeshType>
inline tmp<meshField<scalar,MeshType>> twoPhaseModel::maxRho() const
{
    tmp<meshField<scalar,MeshType>> tMaxRho
    (
        max(this->rho1<MeshType>(), this->rho2<MeshType>())
    );
    tMaxRho->setRestrictionScheme("volumeWeighted");
    return tMaxRho;
}

template<class MeshType>
inline tmp<meshField<scalar,MeshType>> twoPhaseModel::meanRho() const
{
    tmp<meshField<scalar,MeshType>> tMeanRho
    (
        new meshField<scalar,MeshType>
        (
            "meanRho",
            (this->rho1<MeshType>() + this->rho2<MeshType>())/2.0
        )
    );
    tMeanRho->setRestrictionScheme("volumeWeighted");
    return tMeanRho;
}

template<class MeshType>
inline tmp<meshField<lowerFaceScalar,MeshType>> twoPhaseModel::nu() const
{
    tmp<meshField<lowerFaceScalar,MeshType>> tNu
    (
        this->mu<MeshType>()
      / ex::interp(this->rho<MeshType>())
    );
    tNu->setRestrictionScheme("faceAverage");
    return tNu;
}

template<class MeshType>
inline tmp<meshField<scalar,MeshType>> twoPhaseModel::v() const
{
    tmp<meshField<scalar,MeshType>> tV(1.0/this->rho<MeshType>());
    tV->setRestrictionScheme("harmonicVolumeWeighted");
    return tV;
}

}

}

}
