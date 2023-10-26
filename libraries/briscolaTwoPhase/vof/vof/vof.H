#ifndef vof_H
#define vof_H

#include "fvMesh.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const dictionary& dict_;

    // Volume fraction field
    colocatedScalarField alpha_;

private:

    // Disallow default bitwise assignment
    void operator=(const vof&);

public:

    // Void fraction threshold
    static const scalar threshold;

    // Runtime type information
    TypeName("vof");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        vof,
        dictionary,
        (const dictionary& dict, const fvMesh& fvMsh),
        (dict, fvMsh)
    );

    // Construct from dictionary
    vof(const dictionary&, const fvMesh&);

    // Copy constructor
    vof(const vof&);

    // Clone vof
    virtual autoPtr<vof> clone() const = 0;

    // Selector
    static autoPtr<vof> New(const dictionary&, const fvMesh&);

    // Destructor
    virtual ~vof();

    // Check if cast to derived type is possible
    template<class VofType>
    inline bool castable() const
    {
        return dynamic_cast<const VofType*>(this);
    }

    // Cast to derived type, if possible. Returns error otherwise.
    template<class VofType>
    inline VofType& cast()
    {
        VofType* vfPtr = dynamic_cast<VofType*>(this);

        if (!vfPtr)
        {
            FatalErrorInFunction
                << "Vof object cannot be cast to " << VofType::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *vfPtr;
    }

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return access to the volume fraction field
    inline colocatedScalarField& alpha() { return alpha_; }

    // Return const access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Solve the Vof advection equation given a flux field
    virtual void solve(const colocatedFaceScalarField&) = 0;

    // Solve the Vof advection equation given a staggered velocity field
    virtual void solve(const staggeredScalarField&);
};

}

}

}

#endif
