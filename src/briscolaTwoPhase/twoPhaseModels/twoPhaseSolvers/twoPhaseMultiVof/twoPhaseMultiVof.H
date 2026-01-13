#ifndef twoPhaseMultiVof_H
#define twoPhaseMultiVof_H

#include "vof.H"
#include "CCL.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class twoPhaseMultiVof
:
    public BaseModel
{
protected:

    // Multiple marker fields
    PtrList<colocatedScalarField> alphas_;

    // Normal schemes
    PtrList<normalScheme> normalSchemes_;

    // Surface tension schemes
    PtrList<surfaceTensionScheme> surfaceTensionSchemes_;

    // The Vof models
    PtrList<vof> vofs_;

    // Connected Component Labeling
    PtrList<CCL> CCL_;

    // Global CCL field
    colocatedLabelField globalCCL_;

    // Colors field (for visualisation only)
    colocatedLabelField colors_;

    // Global total number of bubbles
    label globalN_;

    // To which alpha field does the tag belong to
    labelList phi_;

    // Connectivity matrix
    List<List<Switch>> connectivityMatrix_;

private:

    // Disallow default bitwise assignment
    void operator=(const twoPhaseMultiVof&);

    // Initialize PtrLists
    void initAlphas();

    // Create new alpha field
    void addAlphaField();

    // Compute global CCL field
    void computeGlobalCCL();

    // Compute global total number of bubbles
    void computeGlobalN();

    // Set phi
    void setPhi();

    // Set global ID's as tags
    void overwriteTagsWithIDs();

    // Set connectivity matrix
    void setConnectivityMatrix();

    // Prevent collisions
    void moveFields();

public:

    // Runtime type information
    TypeName("twoPhaseMultiVof");

    // Construct from IOdictionary
    twoPhaseMultiVof(const fvMesh&, const IOdictionary&);

    // Copy constructor
    twoPhaseMultiVof(const twoPhaseMultiVof&);

    // Clone twoPhaseMultiVof
    virtual autoPtr<twoPhaseModel> clone() const
    {
        return autoPtr<twoPhaseModel>(new twoPhaseMultiVof(*this));
    }

    // Return surface tension flux
    virtual tmp<colocatedScalarFaceField> flux();

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "twoPhaseMultiVof.C"
#endif

#endif
