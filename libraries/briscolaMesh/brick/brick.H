#ifndef brick_H
#define brick_H

#include "meshObject.H"
#include "PtrList.H"
#include "face.H"
#include "block.H"
#include "geometryFun.H"
#include "dictionary.H"
#include "labelVector.H"
#include "grading.H"
#include "labelVectorBlock.H"
#include "labelTensor.H"

namespace Foam
{

namespace briscola
{

class brick
:
    public meshObject<geometry>
{
    // Reference to the brick dictionary
    const dictionary& dict_;

    // Block of vertices
    labelBlock v_;

    // Number of cells
    labelVector N_;

    // Grading object
    autoPtr<grading> grading_;

    // List of faces
    PtrList<face> faces_;

    // Create or recreate faces
    void createFaces();

public:

    // Construct from components
    brick(const geometry&, const label, const dictionary&);

    // Copy constructor
    brick(const brick&);

    // Clone this brick (needed by PtrList)
    inline autoPtr<brick> clone() const
    {
        return autoPtr<brick>(new brick(*this));
    }

    // Destructor
    virtual ~brick();

    // Access to the geometry
    inline const geometry& parentGeometry() const
    {
        return parent_;
    }

    // Access to the dictionary
    inline const dictionary& dict() const
    {
        return dict_;
    }

    // Number of cells
    inline labelVector N() const
    {
        return N_;
    }

    // Access to the list of faces
    inline const PtrList<face>& faces() const
    {
        return faces_;
    }

    // Access to the brick transformation tensor
    inline labelTensor& T() { return v_.T(); }
    inline const labelTensor& T() const { return v_.T(); }

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

    // Zeroth face (x-)
    inline const face& f0() const
    {
        return faces_[0];
    }

    // First face (x+)
    inline const face& f1() const
    {
        return faces_[1];
    }

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

    // Second face (y-)
    inline const face& f2() const
    {
        return faces_[2];
    }

    // Third face (y+)
    inline const face& f3() const
    {
        return faces_[3];
    }

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

    // Fourth face (z-)
    inline const face& f4() const
    {
        return faces_[4];
    }

    // Fifth face (z+)
    inline const face& f5() const
    {
        return faces_[5];
    }

    // Nth face
    inline const face& f(const label i) const
    {
        return faces_[i];
    }

    //  Edges in x-direction numbering
    //
    //         # # # # (1) # # #
    //         # #             . #
    //         #   #           .   #
    //         #     #         .     #
    //         #       # # # # (3) # # #
    //         # . . . # (0) . .       #
    //           #     #         .     #
    //             #   #           .   #
    //               # #             . #
    //                 # # # # (2) # # #

    // Zeroth edge (x- to x+ for y-,z-)
    inline const edge& e0() const
    {
        return f4().e0();
    }

    // First edge (x- to x+ for y+,z-)
    inline const edge& e1() const
    {
        return f4().e1();
    }

    // Second edge (x- to x+ for y-,z+)
    inline const edge& e2() const
    {
        return f5().e0();
    }

    // Third edge (x- to x+ for y+,z+)
    inline const edge& e3() const
    {
        return f5().e1();
    }

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
    //

    // Fourth edge (y- to y+ for x-,z-)
    inline const edge& e4() const
    {
        return f4().e2();
    }

    // Fifth edge (y- to y+ for x+,z-)
    inline const edge& e5() const
    {
        return f4().e3();
    }

    // Sixth edge (y- to y+ for x-,z+)
    inline const edge& e6() const
    {
        return f5().e2();
    }

    // Seventh edge (y- to y+ for x+,z+)
    inline const edge& e7() const
    {
        return f5().e3();
    }

    //  Edges in z-direction numbering
    //
    //  # # # # # # # # #
    //  # #             . #
    //  #  (10)         .  (11)
    //  #     #         .     #
    //  #       # # # # # # # # #
    //  # . . . # . . . .       #
    //    #     #         .     #
    //     (8)  #          (9)  #
    //        # #             . #
    //          # # # # # # # # #
    //

    // Eighth edge (z- to z+ for x-,y-)
    inline const edge& e8() const
    {
        return f2().e2();
    }

    // Nineth edge (z- to z+ for x+,y-)
    inline const edge& e9() const
    {
        return f2().e3();
    }

    // Tenth edge (z- to z+ for x-,y+)
    inline const edge& e10() const
    {
        return f3().e2();
    }

    // Eleventh edge (z- to z+ for x+,y+)
    inline const edge& e11() const
    {
        return f3().e3();
    }

    // Nth edge
    inline const edge& e(const label i) const
    {
        return
            i == 0  ? e0()
          : i == 1  ? e1()
          : i == 2  ? e2()
          : i == 3  ? e3()
          : i == 4  ? e4()
          : i == 5  ? e5()
          : i == 6  ? e6()
          : i == 7  ? e7()
          : i == 8  ? e8()
          : i == 9  ? e9()
          : i == 10 ? e10()
          :           e11();
    }

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

    // Zeroth vertex (x-, y-, z-)
    inline const vertex& v0() const
    {
        return e0().v0();
    }

    // First vertex (x+, y-, z-)
    inline const vertex& v1() const
    {
        return e0().v1();
    }

    // Second vertex (x-, y+, z-)
    inline const vertex& v2() const
    {
        return e1().v0();
    }

    // Third vertex (x+, y+, z-)
    inline const vertex& v3() const
    {
        return e1().v1();
    }

    // Fourth vertex (x-, y-, z+)
    inline const vertex& v4() const
    {
        return e2().v0();
    }

    // Fifth vertex (x+, y-, z+)
    inline const vertex& v5() const
    {
        return e2().v1();
    }

    // Sixth vertex (x-, y+, z+)
    inline const vertex& v6() const
    {
        return e3().v0();
    }

    // Seventh vertex (x+, y+, z+)
    inline const vertex& v7() const
    {
        return e3().v1();
    }

    // Nth vertex
    inline const vertex& v(const label i) const
    {
        return
            i == 0 ? v0()
          : i == 1 ? v1()
          : i == 2 ? v2()
          : i == 3 ? v3()
          : i == 4 ? v4()
          : i == 5 ? v5()
          : i == 6 ? v6()
          :          v7();
    }

    // Block of vertices
    inline const labelBlock& v() const
    {
        return v_;
    }

    // Volume of the brick (using tet volumes)
    inline scalar volume() const
    {
        return
            tetVolume(v4(), v5(), v0(), v6())
          + tetVolume(v1(), v0(), v5(), v3())
          + tetVolume(v2(), v3(), v6(), v0())
          + tetVolume(v7(), v6(), v3(), v5())
          + tetVolume(v6(), v0(), v3(), v5());
    }

    // Volume of the brick (using magnitude of tet volumes)
    inline scalar magVolume() const
    {
        return
            Foam::mag(tetVolume(v4(), v5(), v0(), v6()))
          + Foam::mag(tetVolume(v1(), v0(), v5(), v3()))
          + Foam::mag(tetVolume(v2(), v3(), v6(), v0()))
          + Foam::mag(tetVolume(v7(), v6(), v3(), v5()))
          + Foam::mag(tetVolume(v6(), v0(), v3(), v5()));
    }

    // Check if brick is right-handed. The coordinates of the vertices are
    // defined in a right-handed coordinate system. Then we can define three
    // vectors x = v(1,0,0) - v(0,0,0), y = v(0,1,0) - v(0,0,0) and z = v(0,0,1)
    // - v(0,0,0) which satisfy dot(x x y, z) > 0 for a right-handed brick-local
    // coordinate system.
    inline bool rightHanded() const
    {
        const vector x(v1()-v0());
        const vector y(v2()-v0());
        const vector z(v4()-v0());

        return ((x^y) & z) > 0;
    }

    inline bool leftHanded() const
    {
        return !rightHanded();
    }

    // Transform this brick given a transformation tensor
    void transform(const labelTensor);

    // Rotate this brick by the given number of steps around the given axis in
    // counter-clockwise direction and return the transformation tensor
    labelTensor rotate(const label, const label);

    // Transfinite interpolation (TFI) mapping the ungraded unit cube onto
    // graded physical space, given three vectors of coordinates which span a
    // subset of the unit cube
    tmp<vectorBlock> TFI
    (
        const scalarField&,
        const scalarField&,
        const scalarField&
    ) const;

    // Equality operator
    bool operator==(const brick&) const;

    // Unequality operator
    bool operator!=(const brick&) const;

    //- Ostream Operator
    friend Ostream& operator<<(Ostream&, const brick&);
};

}

}

#endif
