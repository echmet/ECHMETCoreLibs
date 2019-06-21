#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_IMPL_HPP
#define ECHMET_CAES_ESTIMATOR_HELPERS_IMPL_HPP

#include "estimator_helpers.h"

namespace ECHMET {
namespace CAES {

template <InstructionSet ISet, bool ThreadSafe>
inline
void calculateDistributionWithDerivative_dbl(const double &v,
					     double *const ECHMET_RESTRICT_PTR distribution,
					     double *const ECHMET_RESTRICT_PTR dDistdV,
					     std::vector<TotalEquilibrium<double, ISet, ThreadSafe>> &totalEquilibria,
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

ECHMET_FORCE_INLINE
void estimateComplexesDistribution_dbl(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				       const size_t totalLigandCopyCount,
				       const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				       double *estimatedConcentrations)
{
	size_t rowCounter = 2;
	size_t ecRowCounter = 2;

	/* Set concentrations of all free forms and ligands first */
	for (const auto cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			estimatedConcentrations[rowCounter] = estConcentrations[ecRowCounter];
			rowCounter++; ecRowCounter++;

			const auto &fv = cn->forms[charge - cn->chargeLow];
			rowCounter += fv.size() - 1;
		}
	}

	memcpy(&estimatedConcentrations[rowCounter], &estConcentrations[ecRowCounter], totalLigandCopyCount * sizeof(double));

	rowCounter += totalLigandCopyCount;
	ecRowCounter += totalLigandCopyCount;

	for (const auto cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const auto &forms = cn->forms[charge - cn->chargeLow];

			for (size_t idx = 1; idx < forms.size(); idx++) {
				const auto f = forms[idx];

				double complexConcentration = X10(f->pB + 3.0) * estimatedConcentrations[f->ancestorGlobalIdx + 2] * estimatedConcentrations[f->ligandIFIdx + LGBlockOffset];
				ECHMET_DEBUG_CODE(fprintf(stderr, "N: %s, myIdx: %zu, GAIdx: %zu, LFIdx: %zu, CC: %g\n", f->name.c_str(), f->myIdx + 2, f->ancestorGlobalIdx + 2, f->ligandIFIdx + LGBlockOffset, CAESRealToDouble(complexConcentration)));
				ECHMET_DEBUG_CODE(fprintf(stderr, "   [GA]=%g, [L]=%g, Kx=%g\n", CAESRealToDouble(estimatedConcentrations[f->ancestorGlobalIdx + 2]),
												 CAESRealToDouble(estimatedConcentrations[f->ligandIFIdx + LGBlockOffset]),
												 CAESRealToDouble(X10(f->pB))));
				estimatedConcentrations[f->myIdx + 2] = complexConcentration;
			}
		}
	}
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_ESTIMATOR_HELPERS_IMPL_HPP
