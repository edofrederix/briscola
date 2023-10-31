#ifndef parabolicFitting_H
#define parabolicFitting_H

#include "curvatureScheme.H"
#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class parabolicFitting
:
    public curvatureScheme
{

    // Boundary type stencil for the LSFIR iteration (true for parallel and
    // periodic boundaries)
    bool boundaryType_[3][3][3];

    // Number of vertex in each face
    int NIPV0_[12] = {0};

    // Vertex index for each face in counter clock wise order
    int IPV0_[12][6] = {0};

    // Number of faces
    int originalNf_;

    // Local volume enforcement problem solver
    LVE lve_;

    // Construct boundaryTypes
    void createBoundaryTypes();

    // Construct NIPV0, IPV0 and originalNf
    void createHexagonDescription();

public:

    // Runtime type information
    TypeName("parabolicFitting");

    // Construct from components
    parabolicFitting
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    parabolicFitting(const parabolicFitting&);

    // Destructor
    virtual ~parabolicFitting();

    // Clone curvatureScheme
    virtual autoPtr<curvatureScheme> clone() const
    {
        return autoPtr<curvatureScheme>(new parabolicFitting(*this));
    }

    // Correct the curvature
    virtual void correct();
};

}

}

}

#endif
