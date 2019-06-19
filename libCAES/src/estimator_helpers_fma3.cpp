#include "estimator_helpers.h"

namespace ECHMET {
namespace CAES {

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, false>> &totalEquilibria,
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
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, true>> &totalEquilibria,
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

} // namespace CAES
} // namespace ECHMET
