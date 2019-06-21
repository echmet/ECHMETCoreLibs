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
	size_t rowCounter = 2;

	for (auto &te : totalEquilibria) {
		double X;
		double dX;
		const auto pack = te.TsAnddTsdV(v, activityCoefficients, X, dX);

		const std::vector<double> & Ts = std::get<0>(pack);
		const std::vector<double> & dTsdV = std::get<1>(pack);

		assert(Ts.size() == dTsdV.size());

		const ECHMETReal c = acRaw[te.concentrationIndex];

		const size_t len = te.len;
		for (size_t idx = 0; idx < len; idx++) {
			const size_t rIdx{rowCounter + idx};
			const double T = Ts[idx];
			const double dT = dTsdV[idx];

			/* Distribution */
			distribution[rIdx] = c * T / X;

			/* dDistdV */
			dDistdV[rIdx] = c * (dT * X - T * dX) / X / X;
		}

		rowCounter += len;
	}
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
	size_t rowCounter = 2;

	for (auto &te : totalEquilibria) {
		double X;
		double dX;
		const auto pack = te.TsAnddTsdV(v, activityCoefficients, X, dX);

		const std::vector<double> & Ts = std::get<0>(pack);
		const std::vector<double> & dTsdV = std::get<1>(pack);

		assert(Ts.size() == dTsdV.size());

		const ECHMETReal c = acRaw[te.concentrationIndex];

		const size_t len = te.len;
		for (size_t idx = 0; idx < len; idx++) {
			const size_t rIdx{rowCounter + idx};
			const double T = Ts[idx];
			const double dT = dTsdV[idx];

			/* Distribution */
			distribution[rIdx] = c * T / X;

			/* dDistdV */
			dDistdV[rIdx] = c * (dT * X - T * dX) / X / X;
		}

		rowCounter += len;
	}
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

} // namespace CAES
} // namespace ECHMET
