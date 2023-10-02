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
    // alpha treshold (takes value from vof)
    scalar threshold_;

    // mesh formed by hexaedrons with planar faces (default true)
    bool planarFaces_;

    // True for using vertex averages instead of centers of mass (default false)
    bool centerAveraging_;

    // Boundary type stencil for LSGIR (true for parallel, periodic and empty boundaries)
    bool boundaryTypeLSGIR_[3][3][3];

    // Boundary type stencil for the LSFIR iteration (true for parallel and periodic boundaries)
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

    // Compute only one component of a vectorial product
    inline scalar vectorialProductProjection
    (
        const vector& v1,
        const vector& v2,
        const int& i
    ) const
    {
        if (i == 0)
        {
            return v1[1] * v2[2] - v1[2] * v2[1];
        }
        else if (i == 1)
        {
            return v1[2] * v2[0] - v1[0] * v2[2];
        }
        else
        {
            return v1[0] * v2[1] - v1[1] * v2[0];
        }
    }

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

    static double time1_;
    static double time2_;
    static double time3_;
    static double time4_;
    static double time5_;
    static double time6_;
    static double time7_;
    static double time8_;
};

}

}

}

#endif
