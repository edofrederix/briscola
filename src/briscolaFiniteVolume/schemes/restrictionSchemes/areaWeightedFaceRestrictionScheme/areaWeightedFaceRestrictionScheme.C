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

    for (int fd = 0; fd < 3; fd++)
        field[fd][l].correctAggData(0);

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
    NotImplemented;
}

}

}

}
