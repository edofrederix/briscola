#include "linearSystemAggregation.H"

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

defineTemplateTypeNameAndDebug(colocatedScalarDiagSystemAggregation, 0);
defineTemplateTypeNameAndDebug(staggeredScalarDiagSystemAggregation, 0);

defineTemplateTypeNameAndDebug(colocatedVectorDiagSystemAggregation, 0);
defineTemplateTypeNameAndDebug(staggeredVectorDiagSystemAggregation, 0);

defineTemplateTypeNameAndDebug(colocatedScalarSystemAggregation, 0);
defineTemplateTypeNameAndDebug(staggeredScalarSystemAggregation, 0);

defineTemplateTypeNameAndDebug(colocatedVectorSystemAggregation, 0);
defineTemplateTypeNameAndDebug(staggeredVectorSystemAggregation, 0);

}

}

}
