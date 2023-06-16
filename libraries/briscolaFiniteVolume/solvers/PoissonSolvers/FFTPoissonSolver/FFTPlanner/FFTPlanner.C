#include "FFTPlanner.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

FFTPlanner::FFTPlanner
(
    const fvMesh& fvMsh
)
:
    I_(fvMsh.msh().decomp().map().legend()[Pstream::nProcs() - 1] + unitXYZ),
    N_(fvMsh.msh().cast<rectilinearMesh>().N())
{
    // Identify initial decomposition type

    decompType_ = 0;

    if (I_.x() > 1 && I_.y() > 1 && I_.z() > 1)
    {
        decompType_ = 1;
    }
    else if (I_.x() == 1 && I_.y() > 1 && I_.z() > 1)
    {
        decompType_ = 2;
    }
    else if (I_.x() > 1 && I_.y() == 1 && I_.z() > 1)
    {
        decompType_ = 3;
    }
    else if (I_.x() > 1 && I_.y() > 1 && I_.z() == 1)
    {
        decompType_ = 4;
    }
    else if (I_.x() == 1 && I_.y() == 1 && I_.z() > 1)
    {
        decompType_ = 5;
    }
    else if (I_.x() == 1 && I_.y() > 1 && I_.z() == 1)
    {
        decompType_ = 6;
    }
    else if (I_.x() > 1 && I_.y() == 1 && I_.z() == 1)
    {
        decompType_ = 7;
    }

    // Select solve direction (for tridiagonal solver)

    const rectilinearMesh& mesh = fvMsh.msh().cast<rectilinearMesh>();

    solveDir_ = -1;

    if (!mesh.uniform().x())
    {
        solveDir_ = 0;
    }
    else if (!mesh.uniform().y())
    {
        solveDir_ = 1;
    }
    else if (!mesh.uniform().z())
    {
        solveDir_ = 2;
    }
    else
    {
        if
        (
            (
                   decompType_ == 0
                || decompType_ == 1
                || decompType_ == 3
                || decompType_ == 4
                || decompType_ == 7
            ) &&
            (
                N_.x() > 1
            )
        )
        {
            solveDir_ = 0;
        }
        else if
        (
            (
                   decompType_ == 2
                || decompType_ == 6
            ) &&
            (
                N_.y() > 1
            )
        )
        {
            solveDir_ = 1;
        }
        else
        {
            solveDir_ = 2;
        }
    }
}

// Destructor

FFTPlanner::~FFTPlanner()
{}

}

}

}
