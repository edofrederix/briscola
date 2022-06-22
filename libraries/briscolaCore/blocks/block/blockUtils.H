#ifndef blockUtils_H
#define blockUtils_H

#include "labelTensor.H"
#include "labelVector.H"
#include "tensor.H"
#include "diagTensor.H"
#include "vector.H"

namespace Foam
{

namespace briscola
{

// Identity

const labelTensor eye(1,0,0,0,1,0,0,0,1);

// Rotationion about x

const labelTensor rotateX1(1,0,0,0,0,-1,0,1,0);
const labelTensor rotateX2(1,0,0,0,-1,0,0,0,-1);
const labelTensor rotateX3(1,0,0,0,0,1,0,-1,0);

const labelTensor rotationsX[] =
{
    eye,
    rotateX1,
    rotateX2,
    rotateX3
};

// Rotationion about y

const labelTensor rotateY1(0,0,1,0,1,0,-1,0,0);
const labelTensor rotateY2(-1,0,0,0,1,0,0,0,-1);
const labelTensor rotateY3(0,0,-1,0,1,0,1,0,0);

const labelTensor rotationsY[] =
{
    eye,
    rotateY1,
    rotateY2,
    rotateY3
};

// Rotationion about z

const labelTensor rotateZ1(0,-1,0,1,0,0,0,0,1);
const labelTensor rotateZ2(-1,0,0,0,-1,0,0,0,1);
const labelTensor rotateZ3(0,1,0,-1,0,0,0,0,1);

const labelTensor rotationsZ[] =
{
    eye,
    rotateZ1,
    rotateZ2,
    rotateZ3
};

// Collection of rotations

const labelTensor rotations[3][4] =
{
    rotationsX,
    rotationsY,
    rotationsZ
};

// One step rotations

const labelTensor rotations1[] =
{
    rotateX1,
    rotateY1,
    rotateZ1
};

// Two step rotations

const labelTensor rotations2[] =
{
    rotateX2,
    rotateY2,
    rotateZ2
};

// Three step rotations

const labelTensor rotations3[] =
{
    rotateX3,
    rotateY3,
    rotateZ3
};

// Reflections

const labelTensor reflectX(-1,0,0,0,1,0,0,0,1);
const labelTensor reflectY(1,0,0,0,-1,0,0,0,1);
const labelTensor reflectZ(1,0,0,0,1,0,0,0,-1);

// Two step rotations

const labelTensor reflections[] =
{
    reflectX,
    reflectY,
    reflectZ
};

// Permutation

const labelTensor permuteXY(0,1,0,1,0,0,0,0,1);
const labelTensor permuteXZ(0,0,1,0,1,0,1,0,0);
const labelTensor permuteYZ(1,0,0,0,0,1,0,1,0);

// Unit label vectors

const labelVector zeroXYZ(0,0,0);

const labelVector unitX(1,0,0);
const labelVector unitY(0,1,0);
const labelVector unitZ(0,0,1);

const labelVector unitXYZ(1,1,1);

const labelVector units[] =
{
    unitX,
    unitY,
    unitZ
};

// Face offsets

const labelVector faceOffset0(-1,0,0);
const labelVector faceOffset1( 1,0,0);

const labelVector faceOffset2(0,-1,0);
const labelVector faceOffset3(0, 1,0);

const labelVector faceOffset4(0,0,-1);
const labelVector faceOffset5(0,0, 1);

const labelVector faceOffsets[] =
{
    faceOffset0,
    faceOffset1,
    faceOffset2,
    faceOffset3,
    faceOffset4,
    faceOffset5
};

// Edge offsets

const labelVector edgeOffset0(0,-1,-1);
const labelVector edgeOffset1(0, 1,-1);
const labelVector edgeOffset2(0,-1, 1);
const labelVector edgeOffset3(0, 1, 1);

const labelVector edgeOffset4(-1,0,-1);
const labelVector edgeOffset5( 1,0,-1);
const labelVector edgeOffset6(-1,0, 1);
const labelVector edgeOffset7( 1,0, 1);

const labelVector edgeOffset8 (-1,-1,0);
const labelVector edgeOffset9 ( 1,-1,0);
const labelVector edgeOffset10(-1, 1,0);
const labelVector edgeOffset11( 1, 1,0);

const labelVector edgeOffsets[] =
{
    edgeOffset0,
    edgeOffset1,
    edgeOffset2,
    edgeOffset3,
    edgeOffset4,
    edgeOffset5,
    edgeOffset6,
    edgeOffset7,
    edgeOffset8,
    edgeOffset9,
    edgeOffset10,
    edgeOffset11
};

// Vertex offsets

const labelVector vertexOffset0(-1,-1,-1);
const labelVector vertexOffset1( 1,-1,-1);
const labelVector vertexOffset2(-1, 1,-1);
const labelVector vertexOffset3( 1, 1,-1);
const labelVector vertexOffset4(-1,-1, 1);
const labelVector vertexOffset5( 1,-1, 1);
const labelVector vertexOffset6(-1, 1, 1);
const labelVector vertexOffset7( 1, 1, 1);

const labelVector vertexOffsets[] =
{
    vertexOffset0,
    vertexOffset1,
    vertexOffset2,
    vertexOffset3,
    vertexOffset4,
    vertexOffset5,
    vertexOffset6,
    vertexOffset7
};

// Global face numbers contained in a block
const label faceNumsInBlock[6] = {0,1,2,3,4,5};

// Global edge numbers contained in faces
const label edgeNumsInFace[6][4] =
{
    {4,6,8,10},
    {5,7,9,11},
    {0,2,8,9},
    {1,3,10,11},
    {0,1,4,5},
    {2,3,6,7}
};

// Global vertex numbers contained in faces
const label vertexNumsInFace[6][4] =
{
    {0,2,4,6},
    {1,3,5,7},
    {0,1,4,5},
    {2,3,6,7},
    {0,1,2,3},
    {4,5,6,7}
};

// Global vertex numbers contained in edges
const label vertexNumsInEdge[12][2] =
{
    {0,1},
    {2,3},
    {4,5},
    {6,7},
    {0,2},
    {1,3},
    {4,6},
    {5,7},
    {0,4},
    {1,5},
    {2,6},
    {3,7}
};

// Return the shift tensor given a rotation tensor
inline labelTensor shift(const labelTensor& R)
{
    labelTensor S(Zero);

    S.xx() = cmptMag(R.x()) != R.x() ? 1 : 0;
    S.yy() = cmptMag(R.y()) != R.y() ? 1 : 0;
    S.zz() = cmptMag(R.z()) != R.z() ? 1 : 0;

    return S;
}

// Default mag is ambiguous for labelVectors
inline scalar mag(const labelVector& a)
{
    return Foam::sqrt(scalar(a.x()*a.x()+a.y()*a.y()+a.z()*a.z()));
}

// Default cmptMultiply is ambiguous for labelVectors
inline labelVector cmptMultiply(const labelVector& a, const labelVector& b)
{
    return labelVector(a.x()*b.x(), a.y()*b.y(), a.z()*b.z());
}

// Default cmptDivide is ambiguous for labelVectors
inline labelVector cmptDivide(const labelVector& a, const labelVector& b)
{
    return labelVector(a.x()/b.x(), a.y()/b.y(), a.z()/b.z());
}

// Return the minimum of each component
inline labelVector cmptMin(const labelVector& a, const labelVector& b)
{
    return labelVector
    (
        Foam::min(a.x(),b.x()),
        Foam::min(a.y(),b.y()),
        Foam::min(a.z(),b.z())
    );
}

// Return the maximum of each component
inline labelVector cmptMax(const labelVector& a, const labelVector& b)
{
    return labelVector
    (
        Foam::max(a.x(),b.x()),
        Foam::max(a.y(),b.y()),
        Foam::max(a.z(),b.z())
    );
}

// Return the face number given a face offset
inline label faceNumber(const labelVector& o)
{
    for (label i = 0; i < 6; i++)
    {
        if (faceOffsets[i] == o)
        {
            return i;
        }
    }

    return -1;
}

// Return the edge number given an edge offset
inline label edgeNumber(const labelVector& o)
{
    for (label i = 0; i < 12; i++)
    {
        if (edgeOffsets[i] == o)
        {
            return i;
        }
    }

    return -1;
}

// Return the vertex number given vertex offset
inline label vertexNumber(const labelVector& o)
{
    for (label i = 0; i < 8; i++)
    {
        if (vertexOffsets[i] == o)
        {
            return i;
        }
    }

    return -1;
}

// Calculate the rotation tensor of one face offset vector to another
inline labelTensor offsetRotation
(
    const labelVector o1,
    const labelVector o2
)
{
    if (mag(o1) != mag(o2))
    {
        FatalErrorInFunction
            << "An offset transformation tensor is only defined for offsets "
            << "with equal length (i.e., face offsets can be transformed to "
            << "face offsets, edge offsets to edge offsets and vertex offsets "
            << "to vertex offsets)." << endl
            << abort(FatalError);
    }

    if (faceNumber(o1) == -1 && edgeNumber(o1) == -1 && vertexNumber(o1) == -1)
    {
        FatalErrorInFunction
            << "Offset vector " << o1 << " is not a face, edge or vertex offset,"
            << " so cannot find offset rotation."
            << endl << abort(FatalError);
    }

    if (faceNumber(o2) == -1 && edgeNumber(o2) == -1 && vertexNumber(o2) == -1)
    {
        FatalErrorInFunction
            << "Offset vector " << o2 << " is not a face, edge or vertex offset,"
            << " so cannot find offset rotation."
            << endl << abort(FatalError);
    }

    // Brute force: just try all combinations of rotations. This could be done
    // much more elegantly, but this approach is the least bug-sensitive, most
    // compact and not too expensive

    for (label i = 0; i < 4; i++)
    for (label j = 0; j < 4; j++)
    for (label k = 0; k < 4; k++)
    {
        const labelTensor T = (rotationsX[i] & rotationsY[j] & rotationsZ[k]);

        if ((T & o1) == o2)
        {
            return T;
        }
    }

    FatalErrorInFunction
        << "Could not find offset rotation tensor." << endl
        << abort(FatalError);

    return eye;
}

// Find the edge number that two faces have in common
inline label sharedEdgeNumber(const label facei, const label facej)
{
    if (facei == facej)
    {
        FatalErrorInFunction
            << "Two equal face numbers given when looking for a shared edge"
            << " number." << endl << abort(FatalError);
    }

    if (facei/2 == facej/2)
    {
        FatalErrorInFunction
            << "Opposing faces have no shared edge."
            << endl << abort(FatalError);
    }

    for (label i = 0; i < 4; i++)
    for (label j = 0; j < 4; j++)
    if (edgeNumsInFace[facei][i] == edgeNumsInFace[facej][j])
    {
        return edgeNumsInFace[facei][i];
    }

    return -1;
}

// Find the vertex number that three faces have in common
inline label sharedVertexNumber
(
    const label facei,
    const label facej,
    const label facek
)
{
    if (facei == facej || facej == facek)
    {
        FatalErrorInFunction
            << "Two equal face numbers given when looking for a shared vertex"
            << " number." << endl << abort(FatalError);
    }

    if (facei/2 == facej/2 || facej/2 == facek/2 || facek/2 == facei/2)
    {
        FatalErrorInFunction
            << "Opposing faces have no shared vertex."
            << endl << abort(FatalError);
    }

    for (label i = 0; i < 4; i++)
    for (label j = 0; j < 4; j++)
    for (label k = 0; k < 4; k++)
    if
    (
        vertexNumsInFace[facei][i] == vertexNumsInFace[facej][j]
     && vertexNumsInFace[facej][j] == vertexNumsInFace[facek][k]
    )
    {
        return vertexNumsInFace[facei][i];
    }

    return -1;
}

}

}

#endif
