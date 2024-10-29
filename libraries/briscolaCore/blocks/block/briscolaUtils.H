#ifndef briscolaUtils_H
#define briscolaUtils_H

#include "labelTensor.H"
#include "labelVector.H"
#include "tensor.H"
#include "diagTensor.H"
#include "vector.H"
#include "scalarList.H"

namespace Foam
{

namespace briscola
{

// Identity

const labelTensor eye(1,0,0,0,1,0,0,0,1);

// Rotation about x

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

// Rotation about y

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

// Rotation about z

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

const labelTensor rotations[][4] =
{
    {
        eye,
        rotateX1,
        rotateX2,
        rotateX3
    },
    {
        eye,
        rotateY1,
        rotateY2,
        rotateY3
    },
    {
        eye,
        rotateZ1,
        rotateZ2,
        rotateZ3
    }
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

// Collection of permutations

const labelTensor permutations[] =
{
    permuteXY,
    permuteXZ,
    permuteYZ
};

// Unit label vectors

const labelVector zeroXYZ(0,0,0);

const labelVector unitX(1,0,0);
const labelVector unitY(0,1,0);
const labelVector unitZ(0,0,1);

const labelVector unitXY(1,1,0);
const labelVector unitXZ(1,0,1);
const labelVector unitYZ(0,1,1);

const labelVector unitXYZ(1,1,1);

const labelVector units[] =
{
    unitX,
    unitY,
    unitZ
};

// Face offsets

//  Faces in x-direction numbering
//
//         # # # # # # # # #
//         # #             . #
//         #   #           .   #
//         #     #         .     #
//  (0) -------o   # # # # # # # # #
//         # . . . # . . . .   o---#----- (1)
//           #     #         .     #
//             #   #           .   #
//               # #             . #
//                 # # # # # # # # #

const labelVector faceOffset0(-1,0,0);
const labelVector faceOffset1( 1,0,0);

//  Faces in y-direction numbering
//
//             (3)
//              |
//  # # # # # # | # #
//  # #         |   . #
//  #   #       o   .   #
//  #     #         .     #
//  #       # # # # # # # # #
//  # . . . # . . . .       #
//    #     #         .     #
//      #   #     o     .   #
//        # #     |       . #
//          # # # # # # # # #
//                |
//               (2)

const labelVector faceOffset2(0,-1,0);
const labelVector faceOffset3(0, 1,0);

//  Faces in z-direction numbering
//
//      (4)
//        \.
//  # # # # # # # # #
//  # #     \       . #
//  #   #    o      .   #
//  #     #         .     #
//  #       # # # # # # # # #
//  # . . . # . . . .       #
//    #     #      o  .     #
//      #   #       \   .   #
//        # #        \    . #
//          # # # # # \ # # #
//                     \.
//                      (5)

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

const labelVector lowerFaceOffsets[] =
{
    faceOffset0,
    faceOffset2,
    faceOffset4
};

const labelVector upperFaceOffsets[] =
{
    faceOffset1,
    faceOffset3,
    faceOffset5
};

// Edge offsets

//  Edges in x-direction numbering
//
//   # # # # (1) # # #
//   # #             . #
//   #   #           .   #
//   #     #         .     #
//   #       # # # # (3) # # #
//   # . . . # (0) . .       #
//     #     #         .     #
//       #   #           .   #
//         # #             . #
//           # # # # (2) # # #

const labelVector edgeOffset0(0,-1,-1);
const labelVector edgeOffset1(0, 1,-1);
const labelVector edgeOffset2(0,-1, 1);
const labelVector edgeOffset3(0, 1, 1);

//  Edges in y-direction numbering
//
//   # # # # # # # # #
//   # #             . #
//   #   #           .   #
//  (4)    #        (5)    #
//   #       # # # # # # # # #
//   # . . . # . . . .       #
//     #     #         .     #
//       #  (6)          .  (7)
//         # #             . #
//           # # # # # # # # #

const labelVector edgeOffset4(-1,0,-1);
const labelVector edgeOffset5( 1,0,-1);
const labelVector edgeOffset6(-1,0, 1);
const labelVector edgeOffset7( 1,0, 1);

//  Edges in z-direction numbering
//
//   # # # # # # # # #
//   # #             . #
//   #  (10)         .  (11)
//   #     #         .     #
//   #       # # # # # # # # #
//   # . . . # . . . .       #
//     #     #         .     #
//      (8)  #          (9)  #
//         # #             . #
//           # # # # # # # # #

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

//  Vertex numbering:
//
//    (2)                 (3)
//      \                 /
//       o # # # # # # # o
//       # #             . #
//       #   #    (6)    .   #    (7)
//       #     #  /      .     #  /
//       #       o # # # # # # # o
//  (0)--o . . . # . . . o       #
//         #     #      /  .     #
//           #   #    (1)    .   #
//             # #             . #
//               o # # # # # # # o
//              /                 \.
//            (4)                 (5)

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

// Unit cube
const vector unitCube[8] =
{
    vector(0,0,0),
    vector(1,0,0),
    vector(0,1,0),
    vector(1,1,0),
    vector(0,0,1),
    vector(1,0,1),
    vector(0,1,1),
    vector(1,1,1)
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

// Global edge numbers contained in faces in counter-clockwise direction as
// viewed from outside
const label edgeNumsInFaceCC[6][4] =
{
    {8,6,10,4},
    {9,5,11,7},
    {0,9,2,8},
    {3,11,1,10},
    {0,4,1,5},
    {2,7,3,6}
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

// Global vertex numbers contained in faces in counter-clockwise direction as
// viewed from outside
const label vertexNumsInFaceCC[6][4] =
{
    {0,4,6,2},
    {5,1,3,7},
    {0,1,5,4},
    {6,7,3,2},
    {1,0,2,3},
    {4,5,7,6}
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

// Decomposition of a block in tets, used throughout the code mostly for
// geometric calculations. We use six tets such that opposing cell faces are
// split consistently along the same diagonal. Ordering of vertices is such that
// (v0-v3), (v1-v3) and (v2-v3) produce a right-handed base.
const label numberOfTets = 6;
const label tetDecomp[numberOfTets][4]
{
    {4,1,2,0},
    {5,2,6,4},
    {5,1,2,4},
    {2,3,5,6},
    {2,3,1,5},
    {6,3,5,7}
};

// Global vertex numbers contained in general hexahedron faces (decomposing each
// face into two triangles) in counter-clockwise direction as viewed from
// outside. The first two indices correspond to the face diagonal.
const label faceTriangleDecompCC[12][3] =
{
    {4,2,0},
    {2,4,6},
    {3,5,1},
    {5,3,7},
    {1,4,0},
    {4,1,5},
    {6,3,2},
    {3,6,7},
    {2,1,0},
    {1,2,3},
    {5,6,4},
    {6,5,7}
};

// Calculate the linear index from indices and block size

inline label lin
(
    const label& i,
    const label& j,
    const label& k,
    const label& l,
    const label& m,
    const label& n
)
{
    return i*m*n + j*n + k;
}

inline label lin
(
    const label& i,
    const label& j,
    const label& k,
    const labelVector& N
)
{
    return lin(i,j,k,N.x(),N.y(),N.z());
}

inline label lin
(
    const labelVector& ijk,
    const label& l,
    const label& m,
    const label& n
)
{
    return lin(ijk.x(),ijk.y(),ijk.z(),l,m,n);
}

inline label lin
(
    const labelVector& ijk,
    const labelVector& N
)
{
    return lin(ijk.x(), ijk.y(), ijk.z(), N.x(), N.y(), N.z());
}

// Calculate the indices from a linear index and block size

inline labelVector ind
(
    const label& l,
    const label&,
    const label& m,
    const label& n
)
{
    return labelVector
    (
        (l)/(m*n),
        (l/n)%m,
        (l%n)
    );
}

inline labelVector ind
(
    const label& l,
    const labelVector& N
)
{
    return ind(l, N.x(), N.y(), N.z());
}

// Flip the components of a vector
template<class Type>
Vector<Type> flip(const Vector<Type>& v)
{
    return Vector<Type>(v[2],v[1],v[0]);
}

// Return index of the neighbor across given face

inline labelVector nei
(
    const label& i,
    const label& j,
    const label& k,
    const label& f
)
{
    return labelVector(i,j,k) + faceOffsets[f];
}

inline labelVector nei(const labelVector& ijk, const label& f)
{
    return nei(ijk.x(), ijk.y(), ijk.z(), f);
}

// Return index of lower neighbor across given face direction

inline labelVector lowerNei
(
    const label& i,
    const label& j,
    const label& k,
    const label& fd
)
{
    return labelVector(i,j,k) - units[fd];
}

inline labelVector lowerNei(const labelVector& ijk, const label& fd)
{
    return lowerNei(ijk.x(), ijk.y(), ijk.z(), fd);
}

// Return index of upper neighbor across given face

inline labelVector upperNei
(
    const label& i,
    const label& j,
    const label& k,
    const label& fd
)
{
    return labelVector(i,j,k) + units[fd];
}

inline labelVector upperNei(const labelVector& ijk, const label& fd)
{
    return upperNei(ijk.x(), ijk.y(), ijk.z(), fd);
}

// Check if the given index is odd
inline bool odd(const label i, const label j, const label k)
{
    return (i+j+k) % 2;
}

// Check if the given index is even
inline bool even(const label i, const label j, const label k)
{
    return !odd(i,j,k);
}

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

// Return the face number given a face offset (no checks)
inline label faceNumber(const labelVector& o)
{
    return
        cmptMax(cmptMag(o)) != 1 ? -1
      : cmptSum(cmptMag(o)) != 1 ? -1
      : o.x() != 0 ? (o.x()+1)/2
      : o.y() != 0 ? (o.y()+1)/2 + 2
      :              (o.z()+1)/2 + 4;
}

// Return the edge number given an edge offset (no checks)
inline label edgeNumber(const labelVector& o)
{
    return
        cmptMax(cmptMag(o)) != 1 ? -1
      : cmptSum(cmptMag(o)) != 2 ? -1
      : o.x() == 0 ? (o.y()+1)/2 + (o.z()+1)
      : o.y() == 0 ? (o.x()+1)/2 + (o.z()+1) + 4
      :              (o.x()+1)/2 + (o.y()+1) + 8;
}

// Return the vertex number given vertex offset (no checks)
inline label vertexNumber(const labelVector& o)
{
    return
        cmptMax(cmptMag(o)) != 1 ? -1
      : cmptSum(cmptMag(o)) != 3 ? -1
      : (o.x()+1)/2 + (o.y()+1) + (o.z()+1)*2;
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
            << "Offset vector " << o1 << " is not a face, edge or vertex "
            << "offset, so cannot find offset rotation."
            << endl << abort(FatalError);
    }

    if (faceNumber(o2) == -1 && edgeNumber(o2) == -1 && vertexNumber(o2) == -1)
    {
        FatalErrorInFunction
            << "Offset vector " << o2 << " is not a face, edge or vertex "
            << "offset, so cannot find offset rotation."
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

// Compute only one component of a cross product
inline scalar crossProductComponent
(
    const vector& v1,
    const vector& v2,
    const label i
)
{
    return
        i == 0 ? v1.y()*v2.z() - v1.z()*v2.y()
      : i == 1 ? v1.z()*v2.x() - v1.x()*v2.z()
      :          v1.x()*v2.y() - v1.y()*v2.x();
}

// Volume of a tetrahedron
inline scalar tetVolume
(
    const vector& a,
    const vector& b,
    const vector& c,
    const vector& d
)
{
    return ((a-d) & ((b-d)^(c-d)))/6.0;
}

// Triangle area vector
inline vector triAreaVector
(
    const vector& a,
    const vector& b,
    const vector& c
)
{
    return 0.5*((b-a) ^ (c-a));
}

// Triangle area
inline scalar triArea
(
    const vector& a,
    const vector& b,
    const vector& c
)
{
    return Foam::mag(triAreaVector(a,b,c));
}

// Clean round-off errors. This is useful for accurate comparison of two
// numbers, for example when checking if a point coordinate falls within a cell.
// In some cases the floating point round-off error can cause a point to
// incorrectly fall just outside a cell. This is avoided by rounding the last
// few significant digits.

inline scalar clean(scalar s)
{
    if (s != 0)
    {
        const scalar factor =
            Foam::pow(10.0, 14 - ceil(Foam::log10(Foam::mag(s))));

        s = round(s*factor)/factor;
    }

    return s;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form clean(VectorSpace<Form, Cmpt, Ncmpts> v)
{
    for (direction i = 0; i < Ncmpts; i++)
        v.v_[i] = clean(v.v_[i]);

    return v;
}

// Make a list from a scalar
inline scalarList list(const scalar s)
{
    return scalarList(1,s);
}

// Make a list from a vector
inline scalarList list(const vector& v)
{
    scalarList l(3);
    for (int i = 0; i < 3; i++) l[i] = v[i];
    return l;
}

// Make a type from list

template<class Type>
inline Type listToType(const scalarList& l)
{
    Type s;
    forAll(l, i)
        s[i] = l[i];
    return s;
}

template<>
inline scalar listToType(const scalarList& l)
{
    return l[0];
}

template<>
inline label listToType(const scalarList& l)
{
    return l[0];
}

}

// Return a string representation of a list (defined in the Foam namespace like
// all other name functions)
template<class Type>
word name(const List<Type>& l)
{
    std::ostringstream buf;

    buf << '(';

    if (l.size() > 0)
    {
        buf << Foam::name(l[0]);
        for (direction i = 1; i < l.size(); i++)
            buf << ',' << Foam::name(l[i]);
    }

    buf << ')';

    return buf.str();
}

}

#endif
