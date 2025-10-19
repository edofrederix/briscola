#include "fluxFaceRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void fluxFaceRestrictionScheme<Type,colocated>::restrict
(
    faceField<Type,colocated>& field,
    const label l,
    const label
)
{
    const labelVector R(field.fvMsh()[l+1].R());

    // Sum fluxes from corresponding fine grid faces

    forAllFacesInSpecificDirection(field, fd, l+1, 0, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector fijk(briscola::cmptMultiply(ijk, R));

        Type value = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            value += field[fd](l,0,fijk+o);
        }

        field[fd](l+1,0,ijk) = value;
    }
}

}

}

}
