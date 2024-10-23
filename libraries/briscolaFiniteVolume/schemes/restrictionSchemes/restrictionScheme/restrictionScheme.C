#include "restrictionScheme.H"

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
    const word restrictionSchemeType
)
{
    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(restrictionSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown restriction scheme "
            << restrictionSchemeType << nl << nl
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
    meshDirection<Type,MeshType>& res,
    const tmp<meshDirection<Type,MeshType>>& tf,
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
    meshLevel<Type,MeshType>& res,
    const meshLevel<Type,MeshType>& f,
    const bool scale
)
{
    forAll(res, d)
        restrict(res[d],f[d],scale);
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
    for (label l = 1; l < res.size(); l++)
        restrict(res[l], f[l-1], scale);
}

template<class Type, class MeshType>
void restrictionScheme<Type,MeshType>::restrict
(
    meshField<Type,MeshType>& res,
    const tmp<meshField<Type,MeshType>>& tf,
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
    const bool scale
)
{
    for (label l = 1; l < res.size(); l++)
        restrict(res[l], res[l-1], scale);
}

}

}

}
