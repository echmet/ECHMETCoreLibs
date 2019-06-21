#include "estimator_helpers_impl.hpp"

namespace ECHMET {
namespace CAES {

template <>
void estimateComplexesDistribution<double, double, InstructionSet::SSE2>
				(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				 const LigandVec<double> *const ECHMET_RESTRICT_PTR allLigands,
				 const size_t totalLigandCopyCount,
				 const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				 double *estimatedConcentrations)
{
	(void)allLigands;

	estimateComplexesDistribution_dbl(complexNuclei, totalLigandCopyCount, estConcentrations,
					  LGBlockOffset, estimatedConcentrations);
}

} // namespace CAES
} // namespace ECHMET
