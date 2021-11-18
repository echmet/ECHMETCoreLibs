#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_HPP
#define ECHMET_CAES_ESTIMATOR_HELPERS_HPP

#include "estimator_helpers.h"
#include <xmmintrin.h>

namespace ECHMET {
namespace CAES {

template <typename T>
class FetchType {
public:
	using Type = T&;
	using CType = const T&;
};

template <>
class FetchType<double> {
public:
	using Type = double;
	using CType = const double;
};

template <typename CAESReal>
inline
void calculateActivityCoefficients(const CAESReal &ionicStrength, std::vector<CAESReal> &activityCoefficients, const std::vector<int> &chargesSquared)
{
	const CAESReal isSqrt = VMath::sqrt<CAESReal>(ionicStrength);

	activityCoefficients[0] = 1.0;
	for (size_t charge = 1; charge < activityCoefficients.size(); charge++)
		activityCoefficients[charge] = activityCoefficientInternal(ionicStrength, isSqrt, chargesSquared[charge]);
}

/*!
 * Calculates distribution of concentrations and its derivative of the entire system
 *
 * @param[in] v Variable in the equilibrium equations.
 * @param[in,out] distribution Resulting array of concentrations. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in,out] dDistdV Resulting array of concentration derivatives. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in] totalEquilibria Vector of objects that descibe the given equilibrium.
 * @param[in] acRaw Array of analytical concentrations.
 * @param[in] activityCoefficients Vector of activity coefficients.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void calculateDistributionWithDerivative(const CAESReal &v,
					 CAESReal *const ECHMET_RESTRICT_PTR distribution,
					 CAESReal *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<CAESReal, ISet, ThreadSafe>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<CAESReal> &activityCoefficients)
{
	size_t rowCounter = 2;

	for (auto &te : totalEquilibria) {
		CAESReal X;
		CAESReal dX;
		const auto pack = te.TsAnddTsdV(v, activityCoefficients, X, dX);

		const std::vector<CAESReal> & Ts = std::get<0>(pack);
		const std::vector<CAESReal> & dTsdV = std::get<1>(pack);

		assert(Ts.size() == dTsdV.size());

		const ECHMETReal c = acRaw[te.concentrationIndex];

		const size_t len = te.len;
		for (size_t idx = 0; idx < len; idx++) {
			const size_t rIdx{rowCounter + idx};
			typename FetchType<CAESReal>::CType T = Ts[idx];
			typename FetchType<CAESReal>::CType dT = dTsdV[idx];

			/* Distribution */
			distribution[rIdx] = c * T / X;

			/* dDistdV */
			dDistdV[rIdx] = c * (dT * X - T * dX) / X / X;
		}

		rowCounter += len;
	}
}

/*!
 * Calculates distribution of concentrations of the entire system
 *
 * @param[in] v Variable in the equilibrium equations.
 * @param[in,out] distribution Resulting vector of concentrations. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in] totalEquilibria Vector of objects that descibe the given equilibria.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
inline
void calculateDistribution(const CAESReal &v,
			   CAESReal * ECHMET_RESTRICT_PTR distribution,
			   std::vector<TotalEquilibrium<CAESReal, ISet, ThreadSafe>> &totalEquilibria,
			   const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
			   const std::vector<CAESReal> &activityCoefficients)
{
	distribution += 2;

	for (auto &te : totalEquilibria) {
		CAESReal X;
		const std::vector<CAESReal> & Ts = te.Ts(v, activityCoefficients, X);

		const ECHMETReal c = acRaw[te.concentrationIndex];
		for (const CAESReal &T : Ts) {
			*distribution = c * T / X;
			distribution++;
		}
	}
}

/*!
 * Calculates initial estimate of concentration of all complex forms.
 *
 * @param[in] complexNuclei Pointer to vector of all complex nuclei.
 * @param[in] allLigands Pointer to vector of all ligands.
 * @param[in] totalLigandCopySize Number of ligand ionic forms that can be processed in a single pass
 * @oaram[in] estConcentrations Array of estimated concentrations of all free (=uncomplexed) species.
 * @param[in] LGBlockOffset Offset of the block that contains concentraions of free ligands in the vector of all forms.
 * @param[in,out] estimatedConcentrations Array of ionic concentrations of all species.
 */
template <typename CAESReal, typename OutputReal, InstructionSet ISet>
void estimateComplexesDistribution(const CNVec<CAESReal> *const ECHMET_RESTRICT_PTR complexNuclei,
				   const LigandVec<CAESReal> *const ECHMET_RESTRICT_PTR allLigands,
				   const size_t totalLigandCopyCount,
				   const CAESReal *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				   OutputReal *estimatedConcentrations)
{
	(void)totalLigandCopyCount;
	size_t rowCounter = 2;
	size_t ecRowCounter = 2;

	/* Set concentrations of all free forms and ligands first */
	for (const auto cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			estimatedConcentrations[rowCounter] = CAESRealToECHMETReal<CAESReal, OutputReal>(estConcentrations[ecRowCounter]);
			rowCounter++; ecRowCounter++;

			const auto &fv = cn->forms[charge - cn->chargeLow];
			rowCounter += fv.size() - 1;
		}
	}

	for (const auto l : *allLigands) {
		for (int charge = l->chargeLow; charge <= l->chargeHigh; charge++) {
			estimatedConcentrations[rowCounter] = CAESRealToECHMETReal<CAESReal, OutputReal>(estConcentrations[ecRowCounter]);
			rowCounter++; ecRowCounter++;
		}
	}

	for (const auto cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const auto &forms = cn->forms[charge - cn->chargeLow];

			for (size_t idx = 1; idx < forms.size(); idx++) {
				const auto f = forms[idx];

				CAESReal complexConcentration = X10(f->pB + 3.0) * estimatedConcentrations[f->ancestorGlobalIdx + 2] * estimatedConcentrations[f->ligandIFIdx + LGBlockOffset];
				ECHMET_DEBUG_CODE(fprintf(stderr, "N: %s, myIdx: %zu, GAIdx: %zu, LFIdx: %zu, CC: %g\n", f->name.c_str(), f->myIdx + 2, f->ancestorGlobalIdx + 2, f->ligandIFIdx + LGBlockOffset, CAESRealToDouble(complexConcentration)));
				ECHMET_DEBUG_CODE(fprintf(stderr, "   [GA]=%g, [L]=%g, Kx=%g\n", CAESRealToDouble(estimatedConcentrations[f->ancestorGlobalIdx + 2]),
												 CAESRealToDouble(estimatedConcentrations[f->ligandIFIdx + LGBlockOffset]),
												 CAESRealToDouble(X10(f->pB))));
				estimatedConcentrations[f->myIdx + 2] = CAESRealToECHMETReal<CAESReal, OutputReal>(complexConcentration);
			}
		}
	}
}

template <typename CAESReal, typename OutputReal, InstructionSet ISet>
void setDistributionFast(const CAESReal *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, OutputReal *estimatedConcentrations)
{
	for (size_t idx = 0; idx < count; idx++)
		estimatedConcentrations[idx] = CAESRealToECHMETReal(estConcentrations[idx]);
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_ESTIMATOR_HELPERS_HPP
