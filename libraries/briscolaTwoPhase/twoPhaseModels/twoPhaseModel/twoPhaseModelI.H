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

    return tMinRho;
}

template<class MeshType>
inline tmp<meshField<scalar,MeshType>> twoPhaseModel::maxRho() const
{
    tmp<meshField<scalar,MeshType>> tMaxRho
    (
        max(this->rho1<MeshType>(), this->rho2<MeshType>())
    );

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

    return tNu;
}

}

}

}
