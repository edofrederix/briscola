#include "mappedBoundaryConditionBase.H"
#include "mappedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
mappedBoundaryConditionBase<Type,MeshType>::mappedBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    DirichletBoundaryCondition<Type,MeshType>
    (
        mshField,
        b,
        List<Type>(MeshType::numberOfDirections, Zero)
    )
{
    // Check if boundary supports a mapped boundary condition

    if (!b.castable<mappedBoundary>())
        FatalErrorInFunction
            << "Boundary " << b.name() << " is not a mapped boundary"
            << endl << abort(FatalError);

    // Read averaging options

    if (this->dict().found("setAverages"))
    {
        this->dict().lookup("setAverages") >> setAverages_;
    }
    else
    {
        setAverages_ =
            List<Switch>(MeshType::numberOfDirections, Switch(false));
    }

    forAll(setAverages_, i)
    {
        if (setAverages_[i])
        {
            this->dict().lookup("averages") >> averages_;
        }
    }
}

}

}

}
