#ifndef vof_H
#define vof_H

#include "twoPhaseModel.H"
#include "meshFields.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof
:
    public regIOobject
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const dictionary& dict_;

    // Reference to the normal scheme
    normalScheme& normal_;

    // Reference to the volume fraction field
    colocatedScalarField& alpha_;

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
        (
            const fvMesh& fvMsh,
            const dictionary& dict,
            normalScheme& normal,
            colocatedScalarField& alpha
        ),
        (fvMsh, dict, normal, alpha)
    );

    // Construct from dictionary
    vof(const fvMesh&, const dictionary&, normalScheme&, colocatedScalarField&);

    // Copy constructor
    vof(const vof&);

    // Clone vof
    virtual autoPtr<vof> clone() const = 0;

    // Selector from two-phase model
    static autoPtr<vof> New
    (
        twoPhaseModel&,
        const dictionary&
    );

    // Selector from components
    static autoPtr<vof> New
    (
        const fvMesh&,
        const dictionary&,
        normalScheme&,
        colocatedScalarField&
    );

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

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the normal scheme
    inline normalScheme& normal() { return normal_; }

    // Const access to the normal scheme
    inline const normalScheme& normal() const { return normal_; }

    // Access to the volume fraction field
    inline colocatedScalarField& alpha() { return alpha_; }

    // Const access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Solve the Vof advection equation given a flux field
    virtual void solve(const colocatedLowerFaceScalarField&) = 0;
};

}

}

}

#endif
