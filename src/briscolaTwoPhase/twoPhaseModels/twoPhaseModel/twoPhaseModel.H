#ifndef twoPhaseModel_H
#define twoPhaseModel_H

#include "fvMesh.H"
#include "meshFields.H"
#include "regIOobject.H"
#include "normalScheme.H"
#include "surfaceTensionScheme.H"
#include "vofField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel
:
    public regIOobject
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const IOdictionary& dict_;

    // Volume fraction field
    vofField alpha_;

    // Gravitational acceleration
    vector g_;

    // Surface tension switch to check if surface tension is active
    Switch tension_;

public:

    // Runtime type information
    TypeName("twoPhaseModel");

    // Colocated run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        twoPhaseModel,
        colocated,
        (
            const fvMesh& fvMsh,
            const IOdictionary& dict
        ),
        (fvMsh, dict)
    );

    // Staggered run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        twoPhaseModel,
        staggered,
        (
            const fvMesh& fvMsh,
            const IOdictionary& dict
        ),
        (fvMsh, dict)
    );

    // Construct from IOdictionary
    twoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    twoPhaseModel(const twoPhaseModel&);

    // Clone twoPhaseModel
    virtual autoPtr<twoPhaseModel> clone() const = 0;

    // Selector
    template<class MeshType>
    static autoPtr<twoPhaseModel> New(const fvMesh&, const IOdictionary&);

    // Check if cast to derived model is possible
    template<class Model>
    inline bool castable() const
    {
        return dynamic_cast<const Model*>(this);
    }

    // Cast to derived model, if possible. Returns error otherwise.
    template<class Model>
    inline Model& cast()
    {
        Model* tpmPtr = dynamic_cast<Model*>(this);

        if (!tpmPtr)
        {
            FatalErrorInFunction
                << "Two-phase model cannot be cast to " << Model::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *tpmPtr;
    }

    template<class Model>
    inline const Model& cast() const
    {
        return const_cast<twoPhaseModel&>(*this).cast<Model>();
    }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return const access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Return access to the volume fraction
    inline colocatedScalarField& alpha() { return alpha_; }

    // Return const access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Gravitational acceleration
    const vector& g() const { return g_; }

    // Return the surface tension switch
    Switch tension() const { return tension_; }

    // Return the interfacial momentum flux
    virtual tmp<colocatedScalarFaceField> flux() = 0;

    // Correct the two-phase model
    virtual void correct() = 0;

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }
};

}

}

}

#endif
