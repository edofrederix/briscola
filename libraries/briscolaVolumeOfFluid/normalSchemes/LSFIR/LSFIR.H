#ifndef LSFIR_H
#define LSFIR_H

#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class LSFIR
:
    public normalScheme
{
    // Normal treshold (defaults to 1e-8)
    scalar threshold_;

    // Compute only one component of a vectorial product
    inline scalar vectorialProductProjection
    (
        const vector& v1,
        const vector& v2,
        const int& i
    ) const
    {
        if (i == 0)
        {
            return v1[1] * v2[2] - v1[2] * v2[1];
        }
        else if (i == 1)
        {
            return v1[2] * v2[0] - v1[0] * v2[2];
        }
        else
        {
            return v1[0] * v2[1] - v1[1] * v2[0];
        }
    }

public:

    // Runtime type information
    TypeName("LSFIR");

    // Construct from components
    LSFIR(const vof&, const dictionary&);

    // Copy constructor
    LSFIR(const LSFIR&);

    // Destructor
    virtual ~LSFIR();

    // Clone normalScheme
    virtual autoPtr<normalScheme> clone() const
    {
        return autoPtr<normalScheme>(new LSFIR(*this));
    }

    // Calculate normal
    virtual tmp<colocatedVectorField> operator()();
};

}

}

}

#endif
