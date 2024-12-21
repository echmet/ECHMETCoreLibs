#include "estimator_helpers_impl.hpp"

namespace ECHMET {
namespace CAES {

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX512, false>
		(v, distribution, dDistdV, totalEquilibria, acRaw, activityCoefficients);
}

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX512, false>
		(v, distribution, dDistdV, totalEquilibria, acRaw);
}

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX512, true>
		(v, distribution, dDistdV, totalEquilibria, acRaw, activityCoefficients);
}

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw)
{
	calculateDistributionWithDerivative_dbl<InstructionSet::AVX512, true>
		(v, distribution, dDistdV, totalEquilibria, acRaw);
}

template <>
void estimateComplexesDistribution<double, double, InstructionSet::AVX512>
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
void setDistributionFast<double, double, InstructionSet::AVX512>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations,
			 const size_t count,
			 double *estimatedConcentrations)
{
	setDistributionFast_dbl(estConcentrations, count, estimatedConcentrations);
}

} // namespace CAES
} // namespace ECHMET
