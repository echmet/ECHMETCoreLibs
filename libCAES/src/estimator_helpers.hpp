#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_HPP
#define ECHMET_CAES_ESTIMATOR_HELPERS_HPP

namespace ECHMET {
namespace CAES {

/*!
 * Calculates distribution of concentrations and its derivative of the entire system
 *
 * @param[in] v Variable in the equilibrium equations.
 * @param[in,out] distribution Resulting vector of concentrations. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in,out] dDistdV Resulting vector of concentration derivatives. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in] totalEquilibria Vector of objects that descibe the given equilibria.
 */
template <typename CAESReal, bool ThreadSafe>
void calculateDistributionWithDerivative(const CAESReal &v, SolverVector<CAESReal> &distribution, SolverVector<CAESReal> &dDistdV, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations)
{
	size_t rowCounter = 2;

	for (TotalEquilibriumBase *teb : totalEquilibria) {
		auto *te = static_cast<TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
		CAESReal X = 0.0;
		CAESReal dX = 0.0;
		const std::vector<CAESReal> & Ts = te->Ts(v, X);
		const std::vector<CAESReal> & dTsdV = te->dTsdV(v, dX);

		assert(Ts.size() == dTsdV.size());

		const CAESReal &c = analyticalConcentrations->elem(te->concentrationIndex);

		for (size_t idx = 0; idx < Ts.size(); idx++) {
			const CAESReal &T = Ts[idx];
			const CAESReal &dT = dTsdV[idx];

			/* Distribution */
			const CAESReal fC = T / X;
			distribution(rowCounter) = c * fC;

			/* dDistdV */
			const CAESReal fD = (dT * X - T * dX) / (X * X);
			dDistdV(rowCounter) = c * fD;

			rowCounter++;
		}
	}
}

/*!
 * Calculates distribution of concentrations of the entire system
 *
 * @param[in] v Variable in the equilibrium equations.
 * @param[in,out] distribution Resulting vector of concentrations. The vector has to be resized by the caller to accomodate all individual concentrations.
 * @param[in] totalEquilibria Vector of objects that descibe the given equilibria.
 */
template <typename CAESReal, bool ThreadSafe>
void calculateDistribution(const CAESReal &v, SolverVector<CAESReal> &distribution, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations)
{
	size_t rowCounter = 2;

	for (TotalEquilibriumBase *teb : totalEquilibria) {
		auto *te = static_cast<TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
		CAESReal X = 0.0;
		const std::vector<CAESReal> & Ts = te->Ts(v, X);

		const CAESReal &c = analyticalConcentrations->elem(te->concentrationIndex);
		for (const CAESReal &T : Ts) {
			const CAESReal fC = c * T / X;

			distribution(rowCounter) = fC;

			rowCounter++;
		}
	}
}

template <typename CAESReal, bool ThreadSafe>
CAESReal calcTotalCharge(const SolverVector<CAESReal> &icConcs, const std::vector<TotalEquilibriumBase *> &totalEquilibria)
{
	CAESReal z = 0;
	size_t rowCounter = 2;

	for (const TotalEquilibriumBase *teb : totalEquilibria) {
		const auto *te = static_cast<const TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
		for (int charge = te->numLow; charge <= te->numHigh; charge++)
			z += icConcs(rowCounter++) * charge;

	}

	return z;
}

/*!
 * Calculates initial estimate of concentration of all complex forms.
 *
 * @param[in] complexNuclei Vector of all complex nuclei.
 * @param[in] allLigands Vector of all ligands.
 * @oaram[in] estConcentrations Vector of estimated concentrations of all free (=uncomplexed) species.
 * @param[in] LGBlockOffset Offset of the block that contains concentraions of free ligands in the vector of all forms.
 * @param[in,out] estimatedConcentrations Vector of ionic concentrations of all species.
 */
template <typename CAESReal>
void estimateComplexesDistribution(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands, const SolverVector<CAESReal> &estConcentrations, const size_t LGBlockOffset, SolverVector<CAESReal> &estimatedConcentrations)
{
	size_t rowCounter = 2;
	size_t ecRowCounter = 2;

	/* Set concentrations of all free forms and ligands first */
	for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			estimatedConcentrations(rowCounter) = estConcentrations(ecRowCounter);
			rowCounter++; ecRowCounter++;

			const FormVec<CAESReal> &fv = cn->forms.at(charge - cn->chargeLow);
			rowCounter += fv.size() - 1;
		}
	}
	for (const Ligand<CAESReal> *l : *allLigands) {
		for (int charge = l->chargeLow; charge <= l->chargeHigh; charge++) {
			estimatedConcentrations(rowCounter) = estConcentrations(ecRowCounter);
			rowCounter++; ecRowCounter++;
		}
	}

	for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const FormVec<CAESReal> &forms = cn->forms.at(charge - cn->chargeLow);

			for (size_t idx = 1; idx < forms.size(); idx++) {
				const Form<CAESReal> *f = forms.at(idx);

				CAESReal complexConcentration = X10(f->pB + 3.0) * estimatedConcentrations(f->ancestorGlobalIdx + 2) * estimatedConcentrations(f->ligandIFIdx + LGBlockOffset);
				ECHMET_DEBUG_CODE(fprintf(stderr, "N: %s, myIdx: %zu, GAIdx: %zu, LFIdx: %zu, CC: %g\n", f->name.c_str(), f->myIdx + 2, f->ancestorGlobalIdx + 2, f->ligandIFIdx + LGBlockOffset, CAESRealToDouble(complexConcentration)));
				ECHMET_DEBUG_CODE(fprintf(stderr, "   [GA]=%g, [L]=%g, Kx=%g\n", CAESRealToDouble(estimatedConcentrations(f->ancestorGlobalIdx + 2)),
												 CAESRealToDouble(estimatedConcentrations(f->ligandIFIdx + LGBlockOffset)),
												 CAESRealToDouble(X10(f->pB))));
				estimatedConcentrations(f->myIdx + 2) = complexConcentration;

			}
		}
	}
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_ESTIMATOR_HELPERS_HPP
