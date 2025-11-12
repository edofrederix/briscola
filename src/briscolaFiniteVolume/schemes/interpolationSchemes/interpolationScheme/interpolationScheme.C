#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::interpolationScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::interpolationScheme
(
    const interpolationScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::~interpolationScheme()
{}

template<class Type, class MeshType>
autoPtr<interpolationScheme<Type,MeshType>>
interpolationScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "interpolationSchemes", schemeName);

    word interpolationSchemeType;
    is >> interpolationSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(interpolationSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << interpolationSchemeType
            << nl << nl << "Valid interpolation schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolationScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<Type,staggered>> tInterp =
        meshField<Type,staggered>::New
        (
            "stagInterp("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,staggered>& interp = tInterp.ref();

    forAllCells(interp, d, i, j, k)
        interp(d,i,j,k) =
            0.5*(field(i,j,k) + field(lowerNeighbor(i,j,k,d)));

    return tInterp;
}

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,staggered>> tInterp = stagInterp(tField());

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

template<class Type>
tmp<faceField<Type,staggered>> stagFaceInterp
(
    const meshField<Type,colocated>& field
)
{
    tmp<faceField<Type,staggered>> tInterp =
        faceField<Type,staggered>::New
        (
            "stagFaceInterp("+field.name()+")",
            field.fvMsh()
        );

    faceField<Type,staggered>& interp = tInterp.ref();

    forAllFaces(interp, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        if (d == fd)
        {
            interp[fd](d,ijk) = field(nei);
        }
        else
        {
            interp[fd](d,ijk) =
                0.25
              * (
                    field(ijk)
                  + field(lowerNeighbor(ijk,fd))
                  + field(nei)
                  + field(lowerNeighbor(nei,fd))
                );
        }
    }

    return tInterp;
}

template<class Type>
tmp<faceField<Type,staggered>> stagFaceInterp
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<Type,staggered>> tInterp =
        stagFaceInterp(tField());

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const meshField<Type,staggered>& field
)
{
    tmp<meshField<Type,colocated>> tInterp =
        meshField<Type,colocated>::New
        (
            "coloInterp("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,colocated>& interp = tInterp.ref();

    interp = Zero;

    forAllCells(field, d, i, j, k)
        interp(i,j,k) +=
            0.5/3.0*(field(d,i,j,k) + field(upperNeighbor(i,j,k,d)));

    return tInterp;
}

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,colocated>> tInterp = stagInterp(tField());

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}


template<class Type>
tmp<faceField<Type,colocated>> coloFaceInterp
(
    const meshField<Type,staggered>& field
)
{
    tmp<faceField<Type,colocated>> tInterp =
        faceField<Type,colocated>::New
        (
            "coloFaceInterp("+field.name()+")",
            field.fvMsh()
        );

    faceField<Type,colocated>& interp = tInterp.ref();

    forAllFaces(interp, fd, i, j, k)
        interp[fd](i,j,k) = field(fd,i,j,k);

    return tInterp;
}

template<class Type>
tmp<faceField<Type,colocated>> coloFaceInterp
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<Type,colocated>> tInterp =
        coloFaceInterp(tField());

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

template<class Type>
tmp<faceField<Type,colocated>> coloFaceInterp
(
    const faceField<Type,staggered>& field
)
{
    tmp<faceField<Type,colocated>> tInterp =
        faceField<Type,colocated>::New
        (
            "coloFaceInterp("+field.name()+")",
            field.fvMsh()
        );

    faceField<Type,colocated>& interp = tInterp.ref();

    forAllFaces(interp, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(upperNeighbor(i,j,k,fd));

        interp[fd](ijk) =
            0.5*(field[fd](fd,ijk) + field[fd](fd,nei));
    }

    return tInterp;
}

template<class Type>
tmp<faceField<Type,colocated>> coloFaceInterp
(
    const tmp<faceField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<Type,colocated>> tInterp =
        coloFaceInterp(tField());

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

}

}

}
