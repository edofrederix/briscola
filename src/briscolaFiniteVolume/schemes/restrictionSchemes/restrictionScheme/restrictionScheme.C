#include "restrictionScheme.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
restrictionScheme<Type,MeshType>::restrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
restrictionScheme<Type,MeshType>::restrictionScheme
(
    const restrictionScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
restrictionScheme<Type,MeshType>::~restrictionScheme()
{}

template<class Type, class MeshType>
autoPtr<restrictionScheme<Type,MeshType>>
restrictionScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word type
)
{
    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find
        (
            type == word::null
          ? restrictionScheme<Type,MeshType>::defaultScheme
          : type
        );

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown restriction scheme "  << type << nl << nl
            << "Valid restriction schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<restrictionScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, scheme::nullStream)
    );
}

template<class Type, class MeshType>
void restrictionScheme<Type,MeshType>::restrict
(
    meshLevel<Type,MeshType>& res,
    const tmp<meshLevel<Type,MeshType>>& tf,
    const bool scale
)
{
    restrict(res,tf(),scale);

    if (tf.isTmp())
        tf.clear();
}

template<class Type, class MeshType>
void restrictionScheme<Type,MeshType>::restrict
(
    meshField<Type,MeshType>& res,
    const meshField<Type,MeshType>& f,
    const bool scale
)
{
    for (label l = 0; l < res.size()-1; l++)
        restrict(res[l+1], f[l], scale);
}

template<class Type, class MeshType>
void restrictionScheme<Type,MeshType>::restrict
(
    meshField<Type,MeshType>& res,
    const tmp<meshField<Type,MeshType>>& tf,
    const bool scale
)
{
    if (tf.isTmp())
        const_cast<tmp<meshField<Type,MeshType>>&>(tf)
            ->correctBoundaryConditions();

    restrict(res,tf(),scale);

    if (tf.isTmp())
        tf.clear();
}

template<class Type, class MeshType>
void restrictionScheme<Type,MeshType>::restrict
(
    meshField<Type,MeshType>& res,
    const bool scale
)
{
    for (label l = 0; l < res.size()-1; l++)
        restrict(res[l+1], res[l], scale);
}

}

}

}
