#include "linearSystem.H"

#include "stencil.H"
#include "diagStencil.H"

#include "scalar.H"
#include "vector.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(colocatedScalarDiagSystem, 0);
defineTemplateTypeNameAndDebug(staggeredScalarDiagSystem, 0);

defineTemplateTypeNameAndDebug(colocatedVectorDiagSystem, 0);
defineTemplateTypeNameAndDebug(staggeredVectorDiagSystem, 0);

defineTemplateTypeNameAndDebug(colocatedScalarSystem, 0);
defineTemplateTypeNameAndDebug(staggeredScalarSystem, 0);

defineTemplateTypeNameAndDebug(colocatedVectorSystem, 0);
defineTemplateTypeNameAndDebug(staggeredVectorSystem, 0);

}

}

}
