#ifndef LSFIR_H
#define LSFIR_H

#include "LSGIR.H"
#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class LSFIR
:
    public LSGIR
{
protected:

    // Compute the plane centroid as the center of mass
    vector centroid
    (
        const labelVector&,
        const List<Switch>&,
        const scalar C
    ) const;

    // Compute the plane centroid as vertex average as fallback
    vector centroidFromVertices
    (
        const labelVector&,
        const List<Switch>&,
        const scalar C
    ) const;

public:

    // Runtime type information
    TypeName("LSFIR");

    // Construct from components
    LSFIR(const fvMesh&, const dictionary&, const colocatedScalarField&);

    // Copy constructor
    LSFIR(const LSFIR&);

    // Destructor
    virtual ~LSFIR();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new LSFIR(*this));
    }

    // Correct the normal
    virtual void correct();
};

}

}

}

#endif
