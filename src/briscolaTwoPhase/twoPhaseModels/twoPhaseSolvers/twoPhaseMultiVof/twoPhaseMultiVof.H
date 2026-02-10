#ifndef twoPhaseMultiVof_H
#define twoPhaseMultiVof_H

#include "vofField.H"
#include "SquareMatrix.H"

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
    PtrList<vofField> alphas_;

    // Global tag field
    colocatedLabelField tags_;

    // Colors field (for visualisation only)
    colocatedLabelField colors_;

    // Global total number of particles
    label N_;

    // To which alpha field does the tag belong to
    labelList phi_;

    // Connectivity matrix
    SquareMatrix<Switch> connectivityMatrix_;

private:

    // Initialize PtrLists
    void initAlphas();

    // Create new alpha field
    void addAlphaField();

    // Compute global tagAlgorithm field
    void computeGlobaltagAlgorithm();

    // Compute global total number of particles
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
