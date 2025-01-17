#include "unstructuredMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(unstructuredMesh, 0);
addToRunTimeSelectionTable(mesh, unstructuredMesh, dictionary);

unstructuredMesh::unstructuredMesh(const IOdictionary& dict)
:
    mesh(dict)
{}

unstructuredMesh::unstructuredMesh(autoPtr<mesh>& mshPtr)
:
    mesh(mshPtr(), true)
{
    mshPtr.clear();
}

unstructuredMesh::unstructuredMesh(const unstructuredMesh& msh)
:
    mesh(msh)
{}

unstructuredMesh::unstructuredMesh(unstructuredMesh& msh, bool reuse)
:
    mesh(msh, reuse)
{}

unstructuredMesh::~unstructuredMesh()
{}

}

}
