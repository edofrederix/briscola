#include "structuredMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(structuredMesh, 0);
addToRunTimeSelectionTable(mesh, structuredMesh, dictionary);

void structuredMesh::setGlobalMeshSize()
{
    const brickTopologyMap& map = topology().map();

    List<bool> x(map.l(), false);
    List<bool> y(map.m(), false);
    List<bool> z(map.n(), false);

    N_ = zeroXYZ;

    forAllBlock(map, i, j, k)
    {
        if (!x[i] && map(i,j,k) > -1)
        {
            N_.x() += bricks()[map(i,j,k)].N().x();
            x[i] = true;
        }

        if (!y[j] && map(i,j,k) > -1)
        {
            N_.y() += bricks()[map(i,j,k)].N().y();
            y[j] = true;
        }

        if (!z[k] && map(i,j,k) > -1)
        {
            N_.z() += bricks()[map(i,j,k)].N().z();
            z[k] = true;
        }
    }
}

structuredMesh::structuredMesh(const IOdictionary& dict)
:
    unstructuredMesh(dict)
{
    setGlobalMeshSize();
}

structuredMesh::structuredMesh(autoPtr<mesh>& mshPtr)
:
    unstructuredMesh(mshPtr)
{
    setGlobalMeshSize();
}

structuredMesh::structuredMesh(structuredMesh& msh, bool reuse)
:
    unstructuredMesh(msh, reuse),
    N_(msh.N_)
{}

structuredMesh::~structuredMesh()
{}

}

}
