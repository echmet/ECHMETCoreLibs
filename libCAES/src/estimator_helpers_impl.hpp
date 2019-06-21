#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_IMPL_HPP
#define ECHMET_CAES_ESTIMATOR_HELPERS_IMPL_HPP

#include "estimator_helpers.h"

namespace ECHMET {
namespace CAES {

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
