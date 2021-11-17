#include "estimator_helpers_impl.hpp"

namespace ECHMET {
namespace CAES {

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX, false>
		(v, distribution, dDistdV, totalEquilibria, acRaw, activityCoefficients);
}

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX, true>
		(v, distribution, dDistdV, totalEquilibria, acRaw, activityCoefficients);
}

template <>
void estimateComplexesDistribution<double, double, InstructionSet::AVX>
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

template <>
void setDistributionFast<double, double, InstructionSet::AVX>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations,
			 const size_t count,
			 double *estimatedConcentrations)
{
	setDistributionFast_dbl(estConcentrations, count, estimatedConcentrations);
}


} // namespace CAES
} // namespace ECHMET
