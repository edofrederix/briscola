#include "averageFaceRestrictionScheme.H"

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
void averageFaceRestrictionScheme<Type,colocated>::restrict
(
    faceField<Type,colocated>& field,
    const label l,
    const label
)
{
    const labelVector R(field.fvMsh()[l+1].R());

    // Average of corresponding fine grid faces

    forAllFacesInSpecificDirection(field, fd, l+1, 0, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector fijk(briscola::cmptMultiply(ijk, R));

        label count = 0;
        Type value = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            value += field[fd](l,0,fijk+o);
            count++;
        }

        field[fd](l+1,0,ijk) = value/scalar(count);
    }
}

// Staggered

template<class Type>
void averageFaceRestrictionScheme<Type,staggered>::restrict
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
