#ifndef ECHMET_CAES_CAES_HPP
#define ECHMET_CAES_CAES_HPP

#include "funcs.h"
#include "solvercontextimpl.h"
#include <algorithm>
#include <cassert>

#define ECHMET_IMPORT_INTERNAL
#include <echmetphchconsts.h>

namespace ECHMET {
namespace CAES {

class FastEstimateFailureException : public std::exception {
public:
	using std::exception::exception;
};

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

/*!
 * Calculates the initial estimation of concentration of all species in the system.
 *
 * @param[in] ctx Solver context.
 * @param[in,out] ionicConcentrations Vector of concentrations of all ionic forms.
 *                The vector shall have the expected size.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of the concentration vectors or the \p ctx
 *         pointer is not castable to \p SolverContextImpl .
 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution.
 * @retval RetCode::E_FAST_ESTIMATE_FAILURE Fast estimation failed to find a solution.
 */
template <typename CAESReal>
RetCode estimateDistributionInternal(const CAESReal &cHInitial, SolverContext *ctx, const RealVec *analyticalConcentrations, SolverVector<CAESReal> &estimatedConcentrations,
				     const bool useFastEstimate) noexcept
{
	SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	estimatedConcentrations.resize(ctxImpl->concentrationCount);

	try {
		const SolverVector<CAESReal> estConcentrations = [&]() {
			if (useFastEstimate) {
				if (ctxImpl->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) {
					return estimatepHFast<CAESReal, false>(cHInitial, ctxImpl->totalEquilibria, analyticalConcentrations,
									       ctxImpl->estimatedIonicConcentrations,ctxImpl->dEstimatedIonicConcentrationsdH);
				} else {
					SolverVector<CAESReal> icConcs(ctxImpl->TECount);
					SolverVector<CAESReal> dIcConcsdH(ctxImpl->TECount);
					return estimatepHFast<CAESReal, true>(cHInitial, ctxImpl->totalEquilibria, analyticalConcentrations, icConcs, dIcConcsdH);
				}
			} else {
				if (ctxImpl->options() & SolverContext::Options::DISABLE_THREAD_SAFETY)
					return estimatepHSafe<CAESReal, false>(ctxImpl->totalEquilibria, analyticalConcentrations, ctxImpl->estimatedIonicConcentrations);
				else {
					SolverVector<CAESReal> icConcs(ctxImpl->TECount);
					return estimatepHSafe<CAESReal, true>(ctxImpl->totalEquilibria, analyticalConcentrations, icConcs);
				}
			}
		}();

		ECHMET_DEBUG_CODE(for (int idx = 0; idx < estConcentrations.size(); idx++) {
				  const CAESReal &v = estConcentrations(idx);
				  fprintf(stderr, "estC: %.4g, pX: %.4g\n", CAESRealToDouble(v), CAESRealToDouble(pX(v)));
				  });

		ECHMET_DEBUG_CODE(fprintf(stderr, "Estimated pH = %g\n", CAESRealToDouble(pX(estConcentrations(0)))));
		/* H+ and OH- are expected to be the first and second item in the vector */
		estimatedConcentrations(0) = estConcentrations(0);
		estimatedConcentrations(1) = estConcentrations(1);

		estimateComplexesDistribution<CAESReal>(ctxImpl->complexNuclei, ctxImpl->allLigands, estConcentrations, ctxImpl->allForms->size() + 2, estimatedConcentrations);
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (const FastEstimateFailureException &) {
		return RetCode::E_FAST_ESTIMATE_FAILURE;
	}

	return RetCode::OK;
}

/*
 * Estimates pH of a system by taking only acidobazic equilibria into account.
 * This is a fast variant that uses Newton-Raphson method to achieve fast convergence
 * at the expense of needing a "good enough" initial estimate and the risk of
 * silent failure if the initial estimate is not "good enough".
 *
 * @param[in] Vector of acidobazic equilibria for all considered species
 *
 * @return Vector of concentrations of all ionic forms including \p H+ and \p OH-
 */
template <typename CAESReal, bool ThreadSafe>
SolverVector<CAESReal> estimatepHFast(const CAESReal &cHInitial, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations,
				      SolverVector<CAESReal> &icConcs, SolverVector<CAESReal> &dIcConcsdH)
{
	const CAESReal KW_298 = CAESReal(PhChConsts::KW_298) * 1e6;
	const CAESReal threshold = electroneturalityPrecision<CAESReal>();
	size_t ctr = 0;
	CAESReal cH = cHInitial;
	CAESReal cHNew;

	auto calcTotalCharge = [](const SolverVector<CAESReal> &icConcs, const std::vector<TotalEquilibriumBase *> &totalEquilibria) {
		CAESReal z = 0;
		size_t rowCounter = 2;

		for (const TotalEquilibriumBase *teb : totalEquilibria) {
			const auto *te = static_cast<const TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
			for (int charge = te->numLow; charge <= te->numHigh; charge++)
				z += icConcs(rowCounter++) * charge;

		}

		return z;
	};

	for (;;) {
		calculateDistributionWithDerivative<CAESReal, ThreadSafe>(cH, icConcs, dIcConcsdH, totalEquilibria, analyticalConcentrations);

		CAESReal z = calcTotalCharge(icConcs, totalEquilibria);
		CAESReal dZ = calcTotalCharge(dIcConcsdH, totalEquilibria);

		z += cH - KW_298 / cH;
		dZ += 1.0 + (KW_298 / cH / cH);

		cHNew = cH - z / dZ;

		//fprintf(stderr, "cH %g, z %g, dZ %g, cHNew, %g\n", CAESRealToDouble(cH), CAESRealToDouble(z), CAESRealToDouble(dZ), CAESRealToDouble(cHNew));

		if (cHNew <= 0.0)
			throw FastEstimateFailureException{};


		if (VMath::abs(cHNew - cH) < threshold)
			break;

		cH = cHNew;

		/* Maximum number of iterations exceeded.
		 * It is likely that the solver has run into trouble,
		 * throw an error and let the user deal with it.
		 */
		if (ctr++ > 50)
			throw FastEstimateFailureException{};
	}

	ECHMET_DEBUG_CODE(fprintf(stderr, "cH = %g\n", CAESRealToDouble(cH)));

	icConcs(0) = cHNew;
	icConcs(1) = KW_298 / cHNew;

	return icConcs;
}

/*
 * Estimates pH of a system by taking only acidobazic equilibria into account.
 * This is a safe version that is guaranteed to converge for any \p H<sub>3<sub>O<sup>+</sup>
 * concentration bounded by \p leftWall and \p rightWall parameters.
 *
 * @param[in] Vector of acidobazic equilibria for all considered species
 *
 * @return Vector of concentrations of all ionic forms including \p H+ and \p OH-
 */
template <typename CAESReal, bool ThreadSafe>
SolverVector<CAESReal> estimatepHSafe(std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations, SolverVector<CAESReal> &icConcs)
{
	const CAESReal KW_298 = CAESReal(PhChConsts::KW_298) * 1e6;
	const CAESReal threshold = electroneturalityPrecision<CAESReal>();
	size_t ctr = 0;
	CAESReal cH = 1.0e-4;
	CAESReal leftWall = 0.0;
	CAESReal rightWall = 100000.0;

	auto calcTotalCharge = [](const SolverVector<CAESReal> &icConcs, const std::vector<TotalEquilibriumBase *> &totalEquilibria) {
		CAESReal z = 0;
		size_t rowCounter = 2;

		for (const TotalEquilibriumBase *teb : totalEquilibria) {
			const auto *te = static_cast<const TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
			for (int charge = te->numLow; charge <= te->numHigh; charge++)
				z += icConcs(rowCounter++) * charge;

		}

		return z;
	};

	for (;;) {
		calculateDistribution<CAESReal, ThreadSafe>(cH, icConcs, totalEquilibria, analyticalConcentrations);

		CAESReal z = calcTotalCharge(icConcs, totalEquilibria);

		z += cH - KW_298 / cH;
		//fprintf(stderr, "cH %g, z %g, dZ %g, cHNew, %g\n", CAESRealToDouble(cH), CAESRealToDouble(z), CAESRealToDouble(dZ), CAESRealToDouble(cHNew));

		if (VMath::abs(z) < threshold)
			break;

		/* Maximum number of iterations exceeded, return what we have so far and
		 * hope for the best */
		if (ctr++ > 1000) {
			icConcs(0) = cH;
			icConcs(1) = KW_298 / cH;
			return icConcs;
		}

		if (z > 0)
			rightWall = cH;
		else
			leftWall = cH;

		cH = (rightWall - leftWall) / 2.0 + leftWall;
	}

	ECHMET_DEBUG_CODE(fprintf(stderr, "cH = %g\n", CAESRealToDouble(cH)));

	//fprintf(stderr, "DONE, iters %zu, cH = %g\n", ctr, CAESRealToDouble(cH));

	icConcs(0) = cH;
	icConcs(1) = KW_298 / cH;

	return icConcs;
}

template <typename CAESReal>
RetCode globalDataToInternal(LigandVec<CAESReal> *allLigands, LigandIonicFormVec<CAESReal> *allLigandIFs, CNVec<CAESReal> *cnVec, FormVec<CAESReal> *allForms, const SysComp::IonicFormVec *inGIfVec, const SysComp::ConstituentVec *gcVec) noexcept
{
	if (inGIfVec->size() < 3)	/* Empty system containing just the solvent */
		return RetCode::OK;

	/* Ionic forms vector contans H3O+ and OH - on top.
	 * This would interfere with the rest of the logic so we just
	 * copy the vector without the first two elements... nasty, right? */
	SysComp::IonicFormVec *gIfVec = inGIfVec->slice(2);
	if (gIfVec == nullptr)
		return RetCode::E_NO_MEMORY;

	/* Create a list of ligands first */
	{
		size_t lifCounter = 0;

		try {
			allLigands->reserve(gcVec->size());
		} catch (std::bad_alloc &) {
			gIfVec->destroy();

			return RetCode::E_NO_MEMORY;
		} catch (std::length_error &) {
			gIfVec->destroy();

			return RetCode::E_DATA_TOO_LARGE;
		}

		for (size_t idx = 0; idx < gcVec->size(); idx++) {
			const SysComp::Constituent *ic = gcVec->elem(idx);

			if (ic->ctype == SysComp::ConstituentType::LIGAND) {
				std::vector<CAESReal> pKas;
				Ligand<CAESReal> *l;

				//ECHMET_DEBUG_CODE(fprintf(stderr, "LG %s TC %.4g\n", ic->name->c_str(), CAESRealToDouble(ic->analyticalConcentration->concentration)));

				try {
					pKas.resize(ic->chargeHigh - ic->chargeLow);
					allLigandIFs->reserve(allLigandIFs->size() + (ic->chargeHigh - ic->chargeLow) + 1);
				} catch (std::bad_alloc &) {
					gIfVec->destroy();

					return RetCode::E_NO_MEMORY;
				} catch (std::length_error &) {
					gIfVec->destroy();

					return RetCode::E_DATA_TOO_LARGE;
				}

				for (int charge = ic->chargeLow; charge < ic->chargeHigh; charge++) {
					const size_t iidx = charge - ic->chargeLow;
					pKas[iidx] = (*(ic->pKas))[iidx];
				}

				try {
					l = new Ligand<CAESReal>(std::string(ic->name->c_str()), ic->chargeLow, ic->chargeHigh, ic->analyticalConcentrationIndex, pKas);
				} catch (std::bad_alloc &) {
					gIfVec->destroy();

					return RetCode::E_NO_MEMORY;
				}

				try {
					for (int charge = ic->chargeLow; charge <= ic->chargeHigh; charge++) {
						const std::string name = std::string(ic->name->c_str()) + "(" + std::to_string(charge) + ")";

						LigandIonicForm<CAESReal> *lif = new LigandIonicForm<CAESReal>(name, charge, lifCounter++, l);
						allLigandIFs->emplace_back(lif);
						l->ionicForms.emplace_back(lif);
					}
					allLigands->emplace_back(l);
				} catch (std::bad_alloc &) {
					delete l;
					gIfVec->destroy();

					return RetCode::E_NO_MEMORY;
				}
			}
		}
		allLigands->shrink_to_fit();
	}

	/* Process complex nuclei */
	{
		size_t cnCounter = 0;

		try {
			cnVec->reserve(gcVec->size());
		} catch (std::bad_alloc &) {
			gIfVec->destroy();

			return RetCode::E_NO_MEMORY;
		}

		for (size_t idx = 0; idx < gcVec->size(); idx++) {
			const SysComp::Constituent *ic = gcVec->elem(idx);

			if (ic->ctype == SysComp::ConstituentType::NUCLEUS) {
				ComplexNucleus<CAESReal> *cn;
				std::vector<CAESReal> pKas;

				try {
					pKas.resize(ic->chargeHigh - ic->chargeLow);
				} catch (std::bad_alloc &) {
					gIfVec->destroy();

					return RetCode::E_NO_MEMORY;
				} catch (std::length_error &) {
					gIfVec->destroy();

					return RetCode::E_DATA_TOO_LARGE;
				}
				for (int charge = ic->chargeLow; charge < ic->chargeHigh; charge++) {
					const size_t iidx = charge - ic->chargeLow;
					pKas[iidx] = (*(ic->pKas))[iidx];
				}

				try {
					cn = new ComplexNucleus<CAESReal>(std::string(ic->name->c_str()), ic->chargeLow, ic->chargeHigh, ic->analyticalConcentrationIndex, pKas);
				} catch (std::bad_alloc &) {
					gIfVec->destroy();

					return RetCode::E_NO_MEMORY;
				}

				RetCode tRet;
				try {
					tRet = initializeForms(cn, allForms, ic, gIfVec, allLigandIFs);
					if (tRet != RetCode::OK) {
						delete cn;
						gIfVec->destroy();

						return tRet;
					}
				} catch (std::bad_alloc &) {
					delete cn;
					gIfVec->destroy();

					return tRet;
				}
				cnVec->emplace_back(cn);
				cnCounter++;
			}
		}
		cnVec->shrink_to_fit();
	}
	gIfVec->destroy();

	return RetCode::OK;
}


template <typename CAESReal>
RetCode initializeForms(ComplexNucleus<CAESReal> *cn, FormVec<CAESReal> *allForms, const SysComp::Constituent *ic, const SysComp::IonicFormVec *gIfVec, const LigandIonicFormVec<CAESReal> *allLigandIFs)
{
	/* This is a somewhat tricky part.
	 * We are basically mirroring the global representation of the system into the internal one.
	 * Using the global representation directly in the solver might be possible but since the
	 * global representation does not keep track of its layout internally, we would have to do
	 * a lot of lookups in the solver itself. Converting the global data into internal data structures
	 * allows us to describe the system and its layout in a way we need to process it efficiently */
	size_t fIdx = 0;
	auto getLigandIFIdx = [](const SysComp::Constituent *l, const int charge, const LigandIonicFormVec<CAESReal> *allLigandIFs) -> size_t {
		/* Ligands are identified by name which is the same in the SysComp and internal Ligand class */
		const std::string lName(l->name->c_str());

		for (size_t idx = 0; idx < allLigandIFs->size(); idx++) {
			const LigandIonicForm<CAESReal> *lIF = (*allLigandIFs)[idx];

			if (lIF->base->name == lName && lIF->charge == charge)
				return idx;
		}

		return SIZE_MAX;
	};

	auto getGlobalAncestorIdx = [](const SysComp::IonicForm *iForm, const SysComp::IonicFormVec *gIfVec) -> size_t {
		/* This is based on the assumption that the order of global and internal ionic forms
		 * is the same. */

		const SysComp::IonicForm *ancestor = iForm->ancestor;
		for (size_t idx = 0; idx < gIfVec->size(); idx++) {
			const SysComp::IonicForm *otherIForm = gIfVec->elem(idx);

			if (*(ancestor->name) == *(otherIForm->name))
				return idx;
		}

		return SIZE_MAX;
	};

	auto getAncestorIdx = [](const std::string &name, const std::vector<Form<CAESReal> *> &forms) -> size_t {
		size_t ctr = 0;
		for (const Form<CAESReal> *_f : forms) {
			if (_f->name == name)
				return ctr;

			ctr++;
		}

		return SIZE_MAX;
	};

	auto addContainedLigand = [](ContainedLigandIonicFormVec<CAESReal> &ligandsContained, const ContainedLigandIonicForm<CAESReal> &cl) {
		for (ContainedLigandIonicForm<CAESReal> &icl : ligandsContained) {
			if (icl.lIF->name == cl.lIF->name) {
				icl.count = cl.count;
				return;
			}
		}

		ligandsContained.emplace_back(cl);
	};

	for (int charge = ic->chargeLow; charge <= ic->chargeHigh; charge++) {
		std::vector<Form<CAESReal> *> cnForms;

		for (; fIdx < ic->ionicForms->size(); fIdx++) {
			const SysComp::IonicForm *iForm = ic->ionicForms->elem(fIdx);
			Form<CAESReal> *f;

			/* We have traversed through all ionic forms for the given charge of
			 * the nucleus so lets advance to the next charge stage */
			if (iForm->nucleusCharge != charge)
				break;

			/* This is the prime ionic form upon which all other ionic forms are built */
			if (iForm->containedLigandIFs == nullptr)
				f = new Form<CAESReal>(iForm->name->c_str(), iForm->totalCharge);
			else {
				const size_t ligandIFIdx = getLigandIFIdx(iForm->ligand, iForm->ligandCharge, allLigandIFs);
				const size_t myIdx = allForms->size();
				const size_t ancestorIdx = getAncestorIdx(iForm->ancestor->name->c_str(), cnForms);
				const size_t globalAncestorIdx = getGlobalAncestorIdx(iForm, gIfVec);
				const ContainedLigandIonicForm<CAESReal> cl(iForm->ligandCount, (*allLigandIFs)[ligandIFIdx]);

				if (VMath::isnan(iForm->pB))
					return RetCode::E_MISSING_PB;

				if (ligandIFIdx == SIZE_MAX || ancestorIdx == SIZE_MAX || globalAncestorIdx == SIZE_MAX)
					return RetCode::E_INVALID_COMPOSITION;

				const Form<CAESReal> *ancestor = (*allForms)[globalAncestorIdx];

				ECHMET_DEBUG_CODE(fprintf(stderr, "N:(%s) AN:(%s),  AIDX %zu, GAIDX %zu\n", iForm->name->c_str(), iForm->ancestor->name->c_str(), ancestorIdx, globalAncestorIdx));
				ContainedLigandIonicFormVec<CAESReal> ligandsContained = ancestor->ligandsContained;
				addContainedLigand(ligandsContained, cl);

				f = new Form<CAESReal>(std::string(iForm->name->c_str()), iForm->totalCharge,
						       ligandIFIdx, ligandsContained, ancestorIdx, globalAncestorIdx, myIdx, iForm->pB);
			}
			allForms->emplace_back(f);
			cnForms.emplace_back(f);
		}
		cn->forms.emplace_back(cnForms);
	}

	return RetCode::OK;
}

/*!
 * Precomputes the part of Jacobian that stays constant througout the calculation.
 *
 * @param[in] complexNuclei Vector of all complexNuclei.
 * @param[in] allLigands Vector of all ligands.
 * @param[in] allFormsCount Count of all complex forms and free complex nuclei present in the system.
 * @param[in] allLigandIFsCount Count of all ligand ionic forms present in the system.
 *
 * @return A pointer to \p SolverMatrix<CAESReal> matrix containing the precomputed Jacobian.
 */
template <typename CAESReal>
SolverMatrix<CAESReal> * prepareJacobian(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands,
					 const size_t allFormsCount, const size_t allLigandIFsCount)
{
	const size_t fullDim = allFormsCount + allLigandIFsCount + 2;
	SolverMatrix<CAESReal> *pJx = new SolverMatrix<CAESReal>{fullDim, fullDim};
	SolverMatrix<CAESReal> &Jx = *pJx;

	Jx.setZero(fullDim, fullDim);

	/* Acidobazic equilibira */
	{
		size_t rowCounter = 2;

		/* Complex nuclei */
		for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
			size_t prevFreeCAIdx;

			for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
				const size_t freeCAIdx = rowCounter;
				const size_t chIdx = charge - cn->chargeLow;
				const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);
				const size_t formsCount = chForms.size();

				if (charge == cn->chargeLow)
					rowCounter += formsCount;
				else {
					Jx(rowCounter, prevFreeCAIdx) = -1.0;
					Jx(rowCounter, rowCounter) = 1.0;
					Jx(rowCounter, 0) = -1.0;
					rowCounter += formsCount;
				}

				prevFreeCAIdx = freeCAIdx;
			}
		}

		/* Ligands */
		rowCounter = allFormsCount + 2;
		for (const Ligand<CAESReal> *l : *allLigands) {
			rowCounter++;

			for (int charge = l->chargeLow + 1; charge <= l->chargeHigh; charge++) {
				Jx(rowCounter, rowCounter - 1) = -1.0;
				Jx(rowCounter, rowCounter) = 1.0;
				Jx(rowCounter, 0) = -1.0;
				rowCounter++;
			}
		}
	}

	/* Individial complexation equilibria */
	{
		size_t rowCounter = 2;

		for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
			/* Fork here for all charges */
			for (const std::vector<Form<CAESReal> *> &cnForms : cn->forms) {
				for (const Form<CAESReal> *f : cnForms) {
					if (f->ligandsContained.size() == 0) {
						/* This is the base form so there are no complexation equilibria to solve here */
						rowCounter++;
						continue;
					}

					/* First two columns are used for H+ and OH-, that is why we are adding 2 to everything */
					const size_t ancestorFormIdx = f->ancestorGlobalIdx + 2;
					const size_t complexIdx = f->myIdx + 2;
					const size_t ligandIFIdx = f->ligandIFIdx + allFormsCount + 2;

					Jx(rowCounter, ancestorFormIdx) = -1.0;
					Jx(rowCounter, complexIdx) = 1.0;
					Jx(rowCounter, ligandIFIdx) = -1.0;

					rowCounter++;
				}
			}
		}
	}

	/* Water ionic product */
	Jx(0, 0) = -1.0;
	Jx(0, 1) = -1.0;

	ECHMET_DEBUG_CODE(
		FILE *outfile = fopen("jmat.out", "w");
		if (outfile == nullptr)
			return;

		fprintf(outfile, "H+;");
		fprintf(outfile, "OH-;");
		for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
			for (const std::vector<Form<CAESReal> *> &forms : cn->forms) {
				for (const Form<CAESReal> *f : forms) {
					fprintf(outfile, "%s;", f->name.c_str());
				}
			}
		}
		for (const Ligand<CAESReal> *l : *allLigands) {
			for (const LigandIonicForm<CAESReal> *lIF : l->ionicForms)
				fprintf(outfile, "%s;", lIF->name.c_str());
		}
		fprintf(outfile, "\n");

		for (size_t row = 0; row < fullDim; row++) {
			for (size_t col = 0; col < fullDim; col++) {
				fprintf(outfile, "%.0g;", CAESRealToDouble(Jx(row, col)));
			}
			fprintf(outfile, "\n");
		}

		fclose(outfile);
	)

	return pJx;
}

template <typename CAESReal>
Solver * createSolverInternal(SolverContext *ctx, const NonidealityCorrections corrections) noexcept
{
	SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
		return nullptr;

	return new (std::nothrow) SolverImpl<CAESReal>(ctxImpl, corrections);
}

/*!
 * Creates a solver context for a given system and intializes
 * total and ionic concentrations vectors and their mappings.
 *
 * @param[in,out] ctx Reference to SolverContext to be created by this function.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Not enough memory to create the sovler context.
 * @retval RetCode::E_DATA_TOO_LARGE Amount of data to be processed is too large.
 * @retval RetCode::E_BAD_INPUT Nonsensical input data.
 * @retval RetCode::E_MISSING_PB Complexation constant was not set.
 */
template <typename CAESReal>
RetCode createSolverContextInternal(SolverContext *&ctx, const SolverContext::Options options, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	auto makeTotalEquilibrium = [options](const auto *ct) -> TotalEquilibriumBase * {
		if (options & SolverContext::Options::DISABLE_THREAD_SAFETY)
			return new TotalEquilibrium<CAESReal, false>(ct->chargeLow, ct->chargeHigh, ct->pKas, ct->analyticalConcentrationIndex);
		return new TotalEquilibrium<CAESReal, true>(ct->chargeLow, ct->chargeHigh, ct->pKas, ct->analyticalConcentrationIndex);
	};

	auto countTEForms = [](const auto &TEVec, const auto &caster) {
		int TECount = 0;
		for (const auto *teb : TEVec) {
			const auto *te = caster(teb);
			TECount += te->numHigh - te->numLow + 1;
		}
		return TECount;
	};

	LigandVec<CAESReal> *allLigands = nullptr;
	LigandIonicFormVec<CAESReal> *allLigandIFs = nullptr;
	CNVec<CAESReal> *complexNuclei = nullptr;
	FormVec<CAESReal> *allForms = nullptr;
	SolverMatrix<CAESReal> *preJacobian = nullptr;
	int TECount = 0;
	std::vector<TotalEquilibriumBase *> totalEquilibria;
	RetCode tRet;
	size_t fullDim;

	/* Create internal representation of the system composition */
	try {
		complexNuclei = new CNVec<CAESReal>();
		allForms = new FormVec<CAESReal>();
		allLigands = new LigandVec<CAESReal>();
		allLigandIFs = new LigandIonicFormVec<CAESReal>();
	} catch (std::bad_alloc &) {
		tRet = RetCode::E_NO_MEMORY;

		goto errout_1;
	}

	tRet = globalDataToInternal(allLigands, allLigandIFs, complexNuclei, allForms, chemSystem.ionicForms, chemSystem.constituents);
	if (tRet != RetCode::OK)
		goto errout_2;

	/* Generate all "subjacobians" for each ComplexNucleus */
	try {
		preJacobian = prepareJacobian(complexNuclei, allLigands, allForms->size(), allLigandIFs->size());
	} catch (std::bad_alloc &) {
		tRet = RetCode::E_NO_MEMORY;

		goto errout_3;
	}

	fullDim = allForms->size() + allLigandIFs->size() + 2;

	try {
		totalEquilibria.reserve(complexNuclei->size() + allLigands->size());

		/* Estimate pH of the system by ignoring any complexation */
		for (const ComplexNucleus<CAESReal> *cn : *complexNuclei)
			totalEquilibria.emplace_back(makeTotalEquilibrium(cn));

		for (const Ligand<CAESReal> *l : *allLigands)
			totalEquilibria.emplace_back(makeTotalEquilibrium(l));


		if (options & SolverContext::Options::DISABLE_THREAD_SAFETY)
			TECount = countTEForms(totalEquilibria, [](const TotalEquilibriumBase *teb) { return static_cast<const TotalEquilibrium<CAESReal, false> *>(teb); });
		else
			TECount = countTEForms(totalEquilibria, [](const TotalEquilibriumBase *teb) { return static_cast<const TotalEquilibrium<CAESReal, true> *>(teb); });
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	try {
		ctx = new SolverContextImpl<CAESReal>(allLigands, allLigandIFs, complexNuclei, allForms, preJacobian,
						      fullDim, chemSystem.constituents->size(),
						      std::move(totalEquilibria), options, TECount + 2);
	} catch (std::bad_alloc &) {
		tRet = RetCode::E_NO_MEMORY;

		goto errout_3;
	}

	return RetCode::OK;

errout_3:
	delete preJacobian;

errout_2:
	releasePointerContainer(complexNuclei);
	releasePointerContainer(allForms);
	releasePointerContainer(allLigands);
	releasePointerContainer(allLigandIFs);

	return tRet;

errout_1:
	delete complexNuclei;
	delete allForms;
	delete allLigands;
	delete allLigandIFs;

	return tRet;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CAES_HPP
