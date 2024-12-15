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

template <typename CAESReal>
static
int findAbsoluteMaximumCharge(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands)
{
	int chargeMax = 1;

	/* Nuclei and their complex forms */
	for (const ComplexNucleus<CAESReal> *cn : *complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const size_t chIdx = charge - cn->chargeLow;
			const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);

			if (std::abs(cn->chargeLow) > chargeMax)
				chargeMax = std::abs(cn->chargeLow);

			if (std::abs(cn->chargeHigh) > chargeMax)
				chargeMax = std::abs(cn->chargeHigh);

			for (const Form<CAESReal> *f : chForms) {
				if (std::abs(f->totalCharge) > chargeMax)
					chargeMax = std::abs(f->totalCharge);

			}
		}
	}

	for (const Ligand<CAESReal> *l : *allLigands) {
		if (std::abs(l->chargeLow) > chargeMax)
			chargeMax = std::abs(l->chargeLow);

		if (std::abs(l->chargeHigh) > chargeMax)
			chargeMax = std::abs(l->chargeHigh);
	}

	return chargeMax;
}

template <typename CAESReal>
static
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
static
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
static
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

template <typename CAESReal, InstructionSet ISet>
static
Solver * makeSolverImpl(bool threadUnsafe, SolverContextImpl<CAESReal> *ctxImpl, const Solver::Options options, const NonidealityCorrections corrections)
{
	if (threadUnsafe)
		return new (std::nothrow) SolverImpl<CAESReal, ISet, false>(ctxImpl, options, corrections);
	return new (std::nothrow) SolverImpl<CAESReal, ISet, true>(ctxImpl, options, corrections);
}

template <typename CAESReal>
static
Solver * createSolverInternal(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) noexcept
{
	SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
		return nullptr;

	const bool threadUnsafe = options & Solver::Options::DISABLE_THREAD_SAFETY;

	switch (detectInstructionSet()) {
	case InstructionSet::GENERIC:
		return makeSolverImpl<CAESReal, InstructionSet::GENERIC>(threadUnsafe, ctxImpl, options, corrections);
	case InstructionSet::SSE2:
		return makeSolverImpl<CAESReal, InstructionSet::SSE2>(threadUnsafe, ctxImpl, options, corrections);
	case InstructionSet::AVX:
		return makeSolverImpl<CAESReal, InstructionSet::AVX>(threadUnsafe, ctxImpl, options, corrections);
	case InstructionSet::FMA3:
		return makeSolverImpl<CAESReal, InstructionSet::FMA3>(threadUnsafe, ctxImpl, options, corrections);
	case InstructionSet::AVX512:
#ifdef ECHMET_DISABLE_AVX512
		return makeSolverImpl<CAESReal, InstructionSet::FMA3>(threadUnsafe, ctxImpl, options, corrections);
#else
		return makeSolverImpl<CAESReal, InstructionSet::AVX512>(threadUnsafe, ctxImpl, options, corrections);
#endif // ECHMET_DISABLE_AVX512
	}

	return nullptr;
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
static
RetCode createSolverContextInternal(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	LigandVec<CAESReal> *allLigands = nullptr;
	LigandIonicFormVec<CAESReal> *allLigandIFs = nullptr;
	CNVec<CAESReal> *complexNuclei = nullptr;
	FormVec<CAESReal> *allForms = nullptr;
	SolverMatrix<CAESReal> *preJacobian = nullptr;
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
		ctx = new SolverContextImpl<CAESReal>(allLigands, allLigandIFs, complexNuclei, allForms, preJacobian,
						      fullDim, chemSystem.constituents->size(),
						      findAbsoluteMaximumCharge<CAESReal>(complexNuclei, allLigands));
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
