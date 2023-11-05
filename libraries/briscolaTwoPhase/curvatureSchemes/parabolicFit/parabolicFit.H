#ifndef parabolicFit_H
#define parabolicFit_H

#include "curvatureScheme.H"
#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class parabolicFit
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
    TypeName("parabolicFit");

    // Construct from components
    parabolicFit
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    parabolicFit(const parabolicFit&);

    // Destructor
    virtual ~parabolicFit();

    // Clone curvatureScheme
    virtual autoPtr<curvatureScheme> clone() const
    {
        return autoPtr<curvatureScheme>(new parabolicFit(*this));
    }

    // Correct the curvature
    virtual void correct();
};

}

}

}

#endif
