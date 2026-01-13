#ifndef functionObjects_slottedSphere_H
#define functionObjects_slottedSphere_H

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

class slottedSphere
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
        return vector
        (
          - 2.0*pi*(x.y() - 0.5),
            2.0*pi*(x.x() - 0.5),
            0
        );
    }

    // Check if point is inside shape
    bool inside(const vector x) const
    {
        return
            mag(x - C_) < R_
        && !(
                x.y() < 0.75 - 0.15 + 2.0*0.125
             && x.x() > 0.5 - 0.025
             && x.x() < 0.5 + 0.025
            );
    }

    // Reference volume fraction field
    tmp<colocatedScalarField> alphaRef() const;

public:

    // Runtime type information
    TypeName("slottedSphere");

    // Construct from time and dictionary
    slottedSphere
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