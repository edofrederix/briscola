#ifndef functionObjects_shearedSphere_H
#define functionObjects_shearedSphere_H

#include "briscolaFunctionObject.H"
#include "meshFields.H"
#include "constants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

using ::Foam::mag;
using ::Foam::cos;
using ::Foam::sin;
using ::Foam::sqr;
using ::Foam::sqrt;
using ::Foam::min;
using ::Foam::max;

using ::Foam::constant::mathematical::pi;

class shearedSphere
:
    public briscolaFunctionObject
{
    // Constants
    const scalar T_;
    const scalar R_;
    const vector C_;

    // Simple uniform cell discretization for volume fraction calculation
    labelVector N_;

    // Velocity given a point, time and time step size
    inline vector v(const vector x, const scalar t, const scalar dt) const
    {
        const scalar r = Foam::mag(vector(x.x()-0.5, x.y()-0.5, 0));

        return vector
        (
            cos(pi*(t + 0.5*dt)/T_)
          * sqr(sin(pi*x.x()))
          * sin(2.0*pi*x.y()),

          - cos(pi*(t + 0.5*dt)/T_)
          * sqr(sin(pi*x.y()))
          * sin(2.0*pi*x.x()),

            cos(pi*(t + 0.5*dt)/T_)
          * sqr(1.0 - r/0.5)
        );
    }

    // Check if point is inside shape
    bool inside(const vector x) const
    {
        return mag(x - C_) < R_;
    }

    // Reference volume fraction field
    tmp<colocatedScalarField> alphaRef() const;

public:

    // Runtime type information
    TypeName("shearedSphere");

    // Construct from time and dictionary
    shearedSphere
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end();
};

}

}

}

}

#endif