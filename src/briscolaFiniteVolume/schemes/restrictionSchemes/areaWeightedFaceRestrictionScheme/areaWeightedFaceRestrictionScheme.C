#include "areaWeightedFaceRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void areaWeightedFaceRestrictionScheme<Type,colocated>::restrict
(
    faceField<Type,colocated>& field,
    const label l,
    const label
)
{
    const labelVector R(field.fvMsh()[l+1].R());

    const colocatedScalarFaceField& faf =
        this->fvMsh().template metrics<colocated>().faceAreas();

    // Face area weighted average of corresponding fine grid faces

    forAllFacesInSpecificDirection(field, fd, l+1, 0, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector fijk(briscola::cmptMultiply(ijk, R));

        scalar area = 0.0;
        Type value = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            const scalar a = faf[fd](l,0,fijk+o);

            area += a;
            value += a*field[fd](l,0,fijk+o);
        }

        field[fd](l+1,0,ijk) = value/area;
    }
}

// Staggered

template<class Type>
void areaWeightedFaceRestrictionScheme<Type,staggered>::restrict
(
    faceField<Type,staggered>& field,
    const label l,
    const label d
)
{
    const labelVector R(field.fvMsh()[l+1].R());

    // First, interpolate to colocated faces on the fine level. The fine level
    // colocated x-faces coincide with the coarse level staggered x-faces for
    // all staggering directions, the fine level colocated y-faces with the
    // coarse level staggered y-faces, etc. The fine level colocated x-faces
    // receive from the first staggering direction, the y-faces from the second,
    // etc.

    FastPtrList<meshDirection<Type,colocated>> coloFine(3);
    forAll(coloFine, fd)
        coloFine.set
        (
            fd,
            meshDirection<Type,colocated>::New(field.fvMsh(),l,0)
        );

    forAll(coloFine, fd)
    for (int i = coloFine[fd].I()[0]; i < coloFine[fd].I()[1] + (fd == 0); i++)
    for (int j = coloFine[fd].I()[2]; j < coloFine[fd].I()[3] + (fd == 1); j++)
    for (int k = coloFine[fd].I()[4]; k < coloFine[fd].I()[5] + (fd == 2); k++)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(upperNeighbor(i,j,k,fd));

        coloFine[fd](ijk) =
            0.5*(field[fd](l,fd,ijk) + field[fd](l,fd,nei));
    }

    // Face area weighted average of corresponding colocated fine grid faces

    const colocatedScalarFaceField& faf =
        this->fvMsh().template metrics<colocated>().faceAreas();

    const labelVector S(coloFine[0].I().lower());
    const labelVector E(coloFine[0].I().upper() - unitXYZ);

    forAllFacesInSpecificDirection(field, fd, l+1, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        // The left/bottom/aft colocated fine grid cell contained by the
        // staggered coarse grid cell has an index scaled by R and padding
        // subtracted

        labelVector fijk =
            briscola::cmptMultiply(ijk,R) - staggered::padding[d];

        scalar area = 0.0;
        Type value = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            // Limit the colocated fine grid face index to only use faces that
            // enclose internal cells, effectively applying homogeneous Neumann

            const labelVector upp
            (
                briscola::cmptMin
                (
                    briscola::cmptMax(fijk + o, S),
                    E + units[fd]
                )
            );

            const scalar a = faf[fd](l,0,upp);

            value += a*coloFine[fd](upp);
            area += a;
        }

        field[fd](l+1,d,ijk) = value/area;
    }
}

}

}

}
