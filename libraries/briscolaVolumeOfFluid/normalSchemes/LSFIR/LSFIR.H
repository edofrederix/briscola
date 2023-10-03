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
    // mesh formed by hexaedrons with planar faces (default true)
    bool planarFaces_;

    // True for using vertex averages instead of centers of mass (default false)
    bool centerAveraging_;

    // Boundary type stencil for LSGIR (true for parallel, periodic and empty
    // boundaries)
    bool boundaryTypeLSGIR_[3][3][3];

    // Boundary type stencil for the LSFIR iteration (true for parallel and
    // periodic boundaries)
    bool boundaryTypeLSFIR_[3][3][3];

    // Number of vertex in each face
    int NIPV0_[12] = {0};

    // Vertex index for each face in counter clock wise order
    int IPV0_[12][6] = {0};

    // Number of faces
    int originalNf_;

    // Construct boundaryTypes
    void createBoundaryTypes();

    // Construct NIPV0, IPV0 and originalNf
    void createHexagonDescription();

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
