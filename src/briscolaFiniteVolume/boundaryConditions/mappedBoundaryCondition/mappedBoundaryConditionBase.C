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
        Zero
    ),
    setAverage_(false)
{
    // Check if boundary supports a mapped boundary condition

    if (!b.castable<mappedBoundary>())
        FatalErrorInFunction
            << "Boundary " << b.name() << " is not a mapped boundary"
            << endl << abort(FatalError);

    // Read averaging options

    if (this->dict().found("setAverage"))
        this->dict().lookup("setAverage") >> setAverage_;

    if (setAverage_)
        this->dict().lookup("average") >> average_;
}

}

}

}
