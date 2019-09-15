#ifndef ECHMET_CAES_SOLVERINTERNAL_HPP
#define ECHMET_CAES_SOLVERINTERNAL_HPP

#include "funcs.h"
#include "vecmath/vecmath.h"
#include <internal/phchconsts_calcs.hpp>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined (ECHMET_COMPILER_MINGW) || defined (ECHMET_COMPILER_MSYS)
#include <x86intrin.h>
#else
#include <xmmintrin.h>
#include <immintrin.h>
#endif // ECHMET_COMPILER_

#include <cassert>

namespace ECHMET {
namespace CAES {

/* For debugging purposes only */
#ifdef ECHMET_DEBUG_OUTPUT
#include <iostream>
template <typename CAESReal>
static
std::ostream & operator<<(std::ostream &ostr, const ECHMET::CAES::SolverMatrix<CAESReal> &m)
{
	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			const CAESReal d = m(row, col);
			ostr << d << " ";
		}

		ostr << "\n";
	}

	return ostr;
}
/* END: For debugging purposes only */
#endif // ECHMET_DEBUG_OUTPUT

template <> template <>
void SolverInternal<double, InstructionSet::SSE2>::VectorizedDelogifier<InstructionSet::SSE2>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <> template <>
void SolverInternal<double, InstructionSet::SSE2>::VectorizedLogifier<InstructionSet::SSE2>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <> template <>
void SolverInternal<double, InstructionSet::AVX>::VectorizedDelogifier<InstructionSet::AVX>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <> template <>
void SolverInternal<double, InstructionSet::AVX>::VectorizedLogifier<InstructionSet::AVX>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <> template <>
void SolverInternal<double, InstructionSet::FMA3>::VectorizedDelogifier<InstructionSet::FMA3>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <> template <>
void SolverInternal<double, InstructionSet::FMA3>::VectorizedLogifier<InstructionSet::FMA3>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);

template <typename CAESReal, InstructionSet ISet>
const char * SolverInternal<CAESReal, ISet>::NumericErrorException::infError = "Result of numeric operation is infinity";

template <typename CAESReal, InstructionSet ISet>
const char * SolverInternal<CAESReal, ISet>::NumericErrorException::nanError = "Result of numeric operation is not a number";

template <typename CAESReal, InstructionSet ISet>
SolverInternal<CAESReal, ISet>::NumericErrorException::NumericErrorException(const ExType t) : std::exception(),
	m_t(t)
{
}

template <typename CAESReal, InstructionSet ISet>
const char * SolverInternal<CAESReal, ISet>::NumericErrorException::what() const noexcept
{
	switch (m_t) {
	case ExType::ET_INF:
		return infError;
	case ExType::ET_NAN:
		return nanError;
	default:
		return "Unknown";
	}
}

/*!
 * SolverInternal c-tor
 *
 * @param[in] ctx SolverContext
 */
template <typename CAESReal, InstructionSet ISet>
SolverInternal<CAESReal, ISet>::SolverInternal(const SolverContextImpl<CAESReal> *ctx) :
	NewtonRaphson<CAESReal, ISet>(ctx->preJacobian->rows()),
	m_ctx(ctx),
	m_allForms(ctx->allForms),
	m_allLigandIFs(ctx->allLigandIFs),
	m_allLigands(ctx->allLigands),
	m_complexNuclei(ctx->complexNuclei),
	m_preJacobian(ctx->preJacobian),
	/* This is horrible but Eigen::Map won't let us do this in a sane way */
	m_pCx_raw(AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(ctx->concentrationCount)),
	m_rCx_raw(AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(ctx->concentrationCount)),
	m_pCx(m_pCx_raw, ctx->concentrationCount),
	m_rCx(m_rCx_raw, ctx->concentrationCount),
	m_vecMath(new VecMath<ISet>()),
	m_vecDelog(*m_vecMath, ctx->concentrationCount)
{
	if (m_pCx_raw == nullptr || m_rCx_raw == nullptr) {
		delete m_vecMath;

		AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_pCx_raw);
		AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_rCx_raw);

		throw std::bad_alloc{};
	}

	initializepACoeffs();

	ECHMET_DEBUG_CODE(std::cerr << "Using instruction set " << ISet << "\n");
}

template <typename CAESReal, InstructionSet ISet>
SolverInternal<CAESReal, ISet>::~SolverInternal()
{
	delete m_vecMath;

	AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_pCx_raw);
	AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_rCx_raw);
}

/*!
 * Calculates the values for Fx column vector.
 *
 * @param[in,out] Fx Column vector Fx to be calculated
 * @param[in] pCx Column vector of concentrations calculated by the previous iteration of the computation
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::ACalculateF(typename NewtonRaphson<CAESReal, ISet>::TX &Fx, const typename NewtonRaphson<CAESReal, ISet>::TX &pCx)
{
	const size_t LGBlockOffset = m_allForms->size() + 2;
	const CAESReal pH = CVI(pCx, 0);
	const CAESReal pOH = CVI(pCx, 1);
	size_t rowCounter = 2; /* First two cells contain c(H+) and c(OH-) */

	ECHMET_DEBUG_CODE(fprintf(stderr, "\n*** Iteration #: %zu\n***\n", this->m_iteration));

	validateVector(pCx);

	/* Precalculate actual values of concetrations */
	m_vecDelog(m_rCx_raw, m_pCx_raw);

	/* Zeroize ligand block of Fx */
	for (int idx = LGBlockOffset; idx < Fx.rows(); idx++)
		CVI(Fx, idx) = 0.0;

	CVI(Fx, 1) = 0.0;
	/* Water ionic product */
	CVI(Fx, 0) = 8.0 - pH - pOH;
	/* Ionic strength correction */
	{
		const CAESReal gW = pACoeff(1);

		CVI(Fx, 0) -= 2.0 * gW;
	}

	size_t prevFreeCAIdx;
	for (const ComplexNucleus<CAESReal> *cn : *m_complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const size_t chIdx = charge - cn->chargeLow;
			const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);
			const size_t formsCount = chForms.size();
			const size_t freeCAIdx = rowCounter;

			/* Cell in the Fx vector corresponding to the lowest free ionic form of the
			 * nucleus is used to calculate mass balance for the given nucleus,
			 * cells in the Fx vector for higher free ionic froms are used
			 * to calculate acidobazic equilibira */
			if (charge == cn->chargeLow) { /* Lowest ionic form, calculate mass balance */
				CAESReal resC = (*m_analyticalConcentrations)(cn->analyticalConcentrationIndex);
				size_t iRowCounter = rowCounter;

				//ECHMET_DEBUG_CODE(std::cerr << cc->name << " IIresC=" << resC << "\n");

				/* Iterate over all substances that contain the given nucleus */
				for (int iCharge = cn->chargeLow; iCharge <= cn->chargeHigh; iCharge++) {
					const size_t iFormsCount = cn->forms.at(iCharge - cn->chargeLow).size();

					for (size_t iCtr = 1; iCtr < iFormsCount; iCtr++)
						resC -= CVI(m_rCx, iRowCounter++); /* Subtract concentration of all forms that contain the given nucleus */

					resC -= CVI(m_rCx, iRowCounter++); /* Subtract concentration of free higher ionic form */
					//ECHMET_DEBUG_CODE(std::cerr << cn->name << " --resC=" << resC << "\n");
				}
				ECHMET_DEBUG_CODE(std::cerr << cn->name << " resC=" << resC << "\n");

				CVI(Fx, rowCounter) = resC;
			} else { /* Higher ionic from, calculate acidobazic equilibrium */
				const CAESReal pKa = cn->pKas.at(charge - cn->chargeLow - 1);

				CVI(Fx, rowCounter) = (pKa - 3.0) - pH - CVI(pCx, prevFreeCAIdx) + CVI(pCx, rowCounter);
				{
					const CAESReal gpH = pACoeff(1);
					const CAESReal gALow = pACoeff(charge - 1);
					const CAESReal gAHigh = pACoeff(charge);

					ECHMET_DEBUG_CODE(fprintf(stderr, "gpH = %g, gpALow = %g, gpAHigh = %g, charge %d\n", CAESRealToDouble(gpH), CAESRealToDouble(gALow), CAESRealToDouble(gAHigh), charge));

					CVI(Fx, rowCounter) += -gpH - gALow + gAHigh;
				}
			}
			prevFreeCAIdx = freeCAIdx;

			/* Electroneutrality */
			CVI(Fx, 1) += CVI(m_rCx, rowCounter) * charge;

			rowCounter++;
			/* The following cells are a result of pML(n) - pML(n - 1) - pL - pB = x */
			ECHMET_DEBUG_CODE(std::cerr << "FCnt " << formsCount << "\n")
			int prevComplexCharge = charge; /* Used in ionic strength correction */
			for (size_t ctr = 1; ctr < formsCount; ctr++) {
				const Form<CAESReal> *f = chForms.at(ctr);
				CAESReal caRes = CVI(pCx, rowCounter); /* pML */

				ECHMET_DEBUG_CODE(std::cerr << "[pML=" << rowCounter << "] [FCA=" << freeCAIdx << "] [DAC=" << f->ancestorGlobalIdx + 2 << "] [pL=" << f->ligandIFIdx + LGBlockOffset << "]\n caRes=" << caRes << "\n");

				/* pML(n-1) */
				if (freeCAIdx == f->ancestorGlobalIdx + 2) {
					ECHMET_DEBUG_CODE(std::cerr << "Ancestor from freeCAIdx\n");
					prevComplexCharge = charge;
				}

				caRes -= CVI(pCx, f->ancestorGlobalIdx + 2);

				/* Subtract pL from Fx-value for pML(n) */
				ECHMET_DEBUG_CODE(std::cerr << " pL IDX " << f->ligandIFIdx + LGBlockOffset << " Val = " <<  CVI(pCx, f->ligandIFIdx + LGBlockOffset) << "\n");
				caRes -= CVI(pCx, f->ligandIFIdx + LGBlockOffset);

				//ECHMET_DEBUG_CODE(std::cerr << "Form for CA(" << charge << ")\t" << f << "\n")

				for (const ContainedLigandIonicForm<CAESReal> &cl : f->ligandsContained) {
					const LigandIonicForm<CAESReal> *lif = cl.lIF;
					size_t offset;

					/* Offset in the Fx matrix to calculate mass balance */
					offset = lif->index - (lif->charge - lif->base->chargeLow);

					ECHMET_DEBUG_CODE(std::cerr << lif->name << " -- IF Offset= " << offset + LGBlockOffset << " LGBlockOffset " << LGBlockOffset << "\n");
					ECHMET_DEBUG_CODE(fprintf(stderr, "lIF idx=%lu, lIF charge=%d, lIF chargeBaseLow=%d, lIF count = %d\n", lif->index, lif->charge, lif->base->chargeLow, cl.count));
					//ECHMET_DEBUG_CODE(fprintf(stderr, "PRE: LGMB=%g, CiF=%g\n", CVI(Fx, offset + LGBlockOffset), CVI(m_rCx, rowCounter)));
					CVI(Fx, LGBlockOffset + offset) -= cl.count * CVI(m_rCx, rowCounter);
					//ECHMET_DEBUG_CODE(fprintf(stderr, "POST: LGMB=%g, CiF=%g\n", CVI(Fx, offset + LGBlockOffset), CVI(m_rCx, rowCounter)));
				}
				/* Electroneutrality */
				//ECHMET_DEBUG_CODE(std::cerr << "Form " << f->name << " total charge = " << totalCharge << "\n");
				CVI(Fx, 1) += CVI(m_rCx, rowCounter) * f->totalCharge;

				{
					const int ligandCharge = m_allLigandIFs->at(f->ligandIFIdx)->charge;
					const CAESReal gL = pACoeff(ligandCharge);
					const CAESReal gMLprev = pACoeff(prevComplexCharge);
					const CAESReal gML = pACoeff(f->totalCharge);

					ECHMET_DEBUG_CODE(std::cerr << "gL " << gL << " gMLprev " << gMLprev << " gML " << gML << " | " << ligandCharge << " " << prevComplexCharge << " " << f->totalCharge << "\n");

					caRes += gML - gMLprev - gL;

					prevComplexCharge = charge;
				}

				/* Subtract pB from Fx-value for pML(n) */
				caRes -= f->pB + 3.0;
				CVI(Fx, rowCounter) = caRes;
				rowCounter++;

				ECHMET_DEBUG_CODE(std::cerr << "caRes " << caRes << "\n")
			}

			ECHMET_DEBUG_CODE(std::cerr << "<-- Next CA\n\n")
		}
	}


	/* Complete ligand balance calculations.
	 * Add analytical concentration for the lowest free ionic form as the corresponding
	 * cell in Fx vector contains ligand mass balance.
	 * Calculate acidobazic equilibira for higher free ionic forms */
	rowCounter = 0;
	for (const Ligand<CAESReal> *l : *m_allLigands) {
		for (int charge = l->chargeLow; charge <= l->chargeHigh; charge++) {
			if (charge == l->chargeLow) { /* Lowest form, do mass balance */
				size_t iRowCounter = rowCounter;

				for (int iCharge = l->chargeLow; iCharge <= l->chargeHigh; iCharge++) {
					CVI(Fx, rowCounter + LGBlockOffset) -= CVI(m_rCx, iRowCounter + LGBlockOffset);
					iRowCounter++;
				}

				ECHMET_DEBUG_CODE(fprintf(stderr, "LGCFree = %g, LGcAn = %g, Fx = %g\n", CAESRealToDouble(CVI(m_rCx, rowCounter + LGBlockOffset)),
													 CAESRealToDouble((*m_analyticalConcentrations)(l->analyticalConcentrationIndex)),
													 CAESRealToDouble(CVI(Fx, rowCounter + LGBlockOffset))));
				CVI(Fx, rowCounter + LGBlockOffset) += (*m_analyticalConcentrations)(l->analyticalConcentrationIndex);
			} else { /* Higher form, do acidobazic equilibrium */
				const CAESReal pKa = l->pKas.at(charge - l->chargeLow - 1) - 3.0;

				//ECHMET_DEBUG_CODE(fprintf(stderr, "pH = %f, PREV = %f, CURR = %f, pKa = %f\n", pH, CVI(pCx, rowCounter + LGBlockOffset - 1), CVI(pCx, rowCounter + LGBlockOffset), pKa));

				CVI(Fx, rowCounter + LGBlockOffset) = pKa - pH - CVI(pCx, rowCounter + LGBlockOffset - 1) + CVI(pCx, rowCounter + LGBlockOffset);
				{
					const CAESReal gpH = pACoeff(1);
					const CAESReal gLLow = pACoeff(charge - 1);
					const CAESReal gLHigh = pACoeff(charge);

					ECHMET_DEBUG_CODE(fprintf(stderr, "%s gpH=%f, gLLow=%f, gLHigh=%f, charge=%d\n",l->name.c_str(), CAESRealToDouble(gpH),
																	 CAESRealToDouble(gLLow),
																	 CAESRealToDouble(gLHigh),
																	 charge));

					CVI(Fx, rowCounter + LGBlockOffset) += -gpH - gLLow + gLHigh;
				}
			}

			/* Electroneutrality */
			CVI(Fx, 1) += CVI(m_rCx, rowCounter + LGBlockOffset) * charge;

			rowCounter++;
		}
	}

	/* Complete electroneutrality by adding [H+] and [OH-] */
	CVI(Fx, 1) += CVI(m_rCx, 0) - CVI(m_rCx, 1); //  X10(pH) - X10(pOH);


	ECHMET_DEBUG_CODE(std::cerr << "--- pCx:\n" << pCx << "\n");
	ECHMET_DEBUG_CODE(std::cerr << "--- rCx:\n" << m_rCx << "\n");
	ECHMET_DEBUG_CODE(std::cerr << "--- Fx:\n" << Fx << "\n");

	validateVector(Fx);
}

/*!
 * Completes the Jacobian matrix.
 *
 * @param[in,out] Jx Jacobian matrix to be completed
 * @param[in] pCx Column vector of concentrations calculated by the previous iteration of the computation
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::ACalculateJ(typename NewtonRaphson<CAESReal, ISet>::TM &Jx, const typename NewtonRaphson<CAESReal, ISet>::TX &pCx)
{
	static_cast<void>(pCx);

	assert(Jx.size() == m_preJacobian->size());

	Jx = *m_preJacobian;

	/* NOTE: I am well aware that some of the calculations below can be done in a single for-loop.
	 *       However, the code is pretty cryptic already and squeezing acidobazic equilibira,
	 *       mass balances and electroneutratily into one furry ball of nonsense would make
	 *       the code 100 % write-only. Do not succumb to any temptations to optimize this! */


	/* NOTE 2: ACalculateF is called before ACalculateJ by the NR-Solver. Therefore we can use
	 *         the values from m_rCx matrix as it is recalculated each time in ACalculateF.
	 *         Jacobian generation *will* break should the calling order change! */

	/* Mass balance rows */
	{
		size_t rowCounter = 2;

		/* Complex nuclei */
		for (const ComplexNucleus<CAESReal> *cn : *m_complexNuclei) {
			size_t totalFormsCount = 0;

			for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++)
				totalFormsCount += cn->forms.at(charge - cn->chargeLow).size();

			for (size_t idx = 0; idx < totalFormsCount; idx++)
				Jx(rowCounter, rowCounter + idx) = CVI(m_rCx, rowCounter + idx) * M_LN10;

			rowCounter += totalFormsCount;
		}

		/* Ligands */
		for (size_t idx = 0; idx < m_allLigands->size(); idx++) {
			const Ligand<CAESReal> *l = m_allLigands->at(idx);

			/* Fill the LG-block of the Jacobian */
			size_t ctr = 0;
			for (int charge = l->chargeLow; charge <= l->chargeHigh; charge++) {
				Jx(rowCounter, rowCounter + ctr) = CVI(m_rCx, rowCounter + ctr) * M_LN10;
				ctr++;
			}

			size_t colCounter = 2;
			for (const Form<CAESReal> *f : *m_allForms) {
				//ECHMET_DEBUG_CODE(std::cerr << "Inspecting form " << f->name << "\n");
				for (const LigandIonicForm<CAESReal> *lif : l->ionicForms) {
					ContainedLigandIonicForm<CAESReal> cl; /* TODO: We are copying the ligand form here, passing a const ptr instead would be a bit more efficient */

					if (isLigandIFContained(f->ligandsContained, lif->name, cl)) {
						Jx(rowCounter, colCounter) = cl.count * CVI(m_rCx, colCounter) * M_LN10;
						ECHMET_DEBUG_CODE(std::cerr << "Ligand found: " << cl.lIF->name << " " << rowCounter << " " << colCounter << cl.count << "\n");
					}
				}
				colCounter++;
			}

			rowCounter += ctr;
		}
	}

	/* Electroneutrality */
	{
		size_t colCounter = 2;

		/* Complex nuclei and their complex forms */
		for (const ComplexNucleus<CAESReal> *cn : *m_complexNuclei) {
			for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
				const size_t chIdx = charge - cn->chargeLow;
				const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);

				for (const Form<CAESReal> *f : chForms) {
					Jx(1, colCounter) = -CVI(m_rCx, colCounter) * f->totalCharge * M_LN10;

					colCounter++;
				}
			}
		}

		/* Ligand ionic forms */
		for (const LigandIonicForm<CAESReal> *l : *m_allLigandIFs) {
			Jx(1, colCounter) = -CVI(m_rCx, colCounter) * l->charge * M_LN10;

			colCounter++;
		}

		/* [H+] and [OH-] */
		Jx(1, 0) = -CVI(m_rCx, 0) * M_LN10;
		Jx(1, 1) = CVI(m_rCx, 1) * M_LN10;
	}

	ECHMET_DEBUG_CODE(
		std::cerr << "--- Jx ---\n";
		std::cerr << "H+ OH- ";
		for (const ComplexNucleus<CAESReal> *cn : *m_complexNuclei) {
			for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
				const std::string chName = cn->name + "(" + std::to_string(charge) + ")";
				const size_t chIdx = charge - cn->chargeLow;
				const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);
				const size_t formsCount = chForms.size();

				std::cerr << chName << " ";

				for (size_t ctr = 1; ctr < formsCount; ctr++) {
					const Form<CAESReal> *f = chForms.at(ctr);

					std::cerr << f->name << " ";
				}
			}
		}

		for (const Ligand<CAESReal> *l : *m_allLigands) {
			for (int charge = l->chargeLow; charge <= l->chargeHigh; charge++)
				std::cerr << l->name << "(" << charge << ") ";
		}

		std::cerr << "\n";

		std::cerr << Jx << "\n";
		);
}

/*!
 * Calculates ionic strength of the solution
 *
 * @return Ionic strength of the solution
 */
template<typename CAESReal, InstructionSet ISet>
CAESReal SolverInternal<CAESReal, ISet>::calculateIonicStrength() const
{
	/* NOTE: This function is called after ASolve() of the NR-Solver has updated the current
	 *       concentrations. Matrix m_rCx is recalculated in the process so we can use it to calculate
	 *       ionic strength here. This *will* break should the calling order ever change!
	 */

	CAESReal ionicStrength = 0.0;
	size_t colCounter = 2;

	/* Nuclei and their complex forms */
	for (const ComplexNucleus<CAESReal> *cn : *m_complexNuclei) {
		for (int charge = cn->chargeLow; charge <= cn->chargeHigh; charge++) {
			const size_t chIdx = charge - cn->chargeLow;
			const FormVec<CAESReal> &chForms = cn->forms.at(chIdx);

			for (const Form<CAESReal> *f : chForms) {
				ionicStrength += CVI(m_rCx, colCounter) * m_ctx->chargesSquared[std::abs(f->totalCharge)];

				colCounter++;
			}
		}
	}

	/* Ligands */
	for (const LigandIonicForm<CAESReal> *l : *m_allLigandIFs) {
			ionicStrength += CVI(m_rCx, colCounter) * m_ctx->chargesSquared[std::abs(l->charge)];

		colCounter++;
	}

	/* [H+] and [OH-] */
	ionicStrength += CVI(m_rCx, 0);
	ionicStrength += CVI(m_rCx, 1);

	return 0.5 * ionicStrength / 1000.0;
}

/*!
 * Returns a pointer to \p SolverContext object that stores the system
 * composition data structures.
 */
template <typename CAESReal, InstructionSet ISet>
const SolverContext * SolverInternal<CAESReal, ISet>::context() const noexcept
{
	return m_ctx;
}

/*!
 * Initializes the vector of p-scaled activity coefficients.
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::initializepACoeffs()
{
	m_pACoeffs.resize(m_ctx->chargesSquared.size());

	for (CAESReal &d : m_pACoeffs)
		d = 0.0;
}

/*!
 * Returns total number of iterations needed to solve the system.
 *
 * @return Number of total iterations. Zero if the system has not been
 * solved or no solution has been found.
 */
template <typename CAESReal, InstructionSet ISet>
SolverIterations SolverInternal<CAESReal, ISet>::iterations() const noexcept
{
	SolverIterations iterations;

	iterations.total = m_totalIterations;
	iterations.outer = m_outerIterations;

	return iterations;
}

/*!
 * Returns the p-scaled activity coefficient for a given charge.
 *
 * @param[in] charge Charge to return the p-scaled activity coefficient for.
 *
 * @return p-scaled activity coefficient.
 */
template<typename CAESReal, InstructionSet ISet>
CAESReal SolverInternal<CAESReal, ISet>::pACoeff(const int charge)
{
	return m_pACoeffs[std::abs(charge)];
}

/*!
 * Recaulculates p-scaled activity for the given ionic strength.
 *
 * @param[in] is Ionic strength.
 */
template<typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::recalculatepACoeffs(const CAESReal &is)
{
	/* McInnes approximation */
	const CAESReal sqrtIs = VMath::sqrt(is);

	for (size_t charge = 1; charge < m_pACoeffs.size(); charge++) {
		const int chSq = m_ctx->chargesSquared[charge];

		m_pACoeffs[charge] = pActivityCoefficientInternal(is, sqrtIs, chSq);
	}
}

/*!
 * Solves the system of equations that describe the thermodynamic
 * equilibira of all components of the system.
 *
 * @param[in] analyticalConcentrations Analytical concentrations of constituents in the system.
 * @param[in] estimatedConcentrations Estimated ionic concentrations.
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] iterations Maximum number of iterations to try.
 * @param[in] inIonicStrength Value of ionic strength in <tt>mM/dm<sup>3</sup></tt> to start with.
 *            Ignored unless correction for Debye-Hückel is enabled.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Not enough memory to perform the calculation.
 * @retval RetCode::E_NRS_FAILURE Newton-Rapshon solver encountered an error during calculation.
 * @retval RetCode::E_NRS_NO_CONVERGENCE Newton-Raphson solver failed to converge within the given number of iterations
 * @retval RetCode::E_NRS_STUCK Greatest change of X-value calculated by the Newton-Raphson solver is below the precision threshold.
 * @retval RetCode::E_IS_NO_CONVERGENCE Solver failed to find a solution within the given number of iterations.
 * @retval RetCode::E_INVALID_ARGUMENT Correction for Debye-Hückel was enabled but initial ionic strength was not positive.
 */
template <typename CAESReal, InstructionSet ISet>
RetCode SolverInternal<CAESReal, ISet>::solve(const SolverVector<CAESReal> *analyticalConcentrations, const CAESReal *estimatedConcentrations,
					      const bool isCorrection, const size_t iterations, const CAESReal &inIonicStrength) noexcept
{
	NumericPrecisionSetter<CAESReal> nps{};

	const size_t maxOuterIterations = 100;
	uint32_t totalIterationsCtr = 0;
	uint32_t outerIterationsCtr = 0;

	m_totalIterations = 0;
	m_outerIterations = 0;
	m_finalIonicStrength = 0;
	m_correctForIS = isCorrection;
	m_analyticalConcentrations = analyticalConcentrations;
	CAESReal ionicStrength;
	CAESReal maxChargepACoeff;
	bool ionicStrengthUnstable;

	for (int idx = 0; idx < m_pCx.rows(); idx++)
		CVI(m_pCx, idx) = pX(estimatedConcentrations[idx]);

	if (m_correctForIS) {
		if (inIonicStrength <= 0.0)
			return RetCode::E_INVALID_ARGUMENT;

		ionicStrength = inIonicStrength;
		recalculatepACoeffs(ionicStrength);
	} else
		ionicStrength = 0.0;

	this->maxIterations = iterations;

	do {
		maxChargepACoeff = m_pACoeffs.back();
		try {
			this->ASolve(&m_pCx);
		} catch (NumericErrorException &ex) {
			ECHMET_DEBUG_CODE(fprintf(stderr, "NRS failure: %s\n", ex.what()));
			return RetCode::E_NRS_FAILURE;
		}

		totalIterationsCtr += this->m_iteration;

		switch (this->GetStatus()) {
		case NewtonRaphson<CAESReal, ISet>::Status::SUCCEEDED:
			break;
		case NewtonRaphson<CAESReal, ISet>::Status::NO_CONVERGENCE:
			ECHMET_DEBUG_CODE(fprintf(stderr, "fMax=%g, dxMax=%g\n", CAESRealToDouble(this->GetMaxF()), CAESRealToDouble(this->GetMaxdx())));
			return RetCode::E_NRS_NO_CONVERGENCE;
		case NewtonRaphson<CAESReal, ISet>::Status::STUCK:
			return RetCode::E_NRS_STUCK;
		case NewtonRaphson<CAESReal, ISet>::Status::NO_SOLUTION:
			return RetCode::E_NRS_NO_SOLUTION;
		default:
			return RetCode::E_NRS_FAILURE;
		}

		if (++outerIterationsCtr > maxOuterIterations)
			return RetCode::E_IS_NO_CONVERGENCE;

		ionicStrength = calculateIonicStrength();
		if (m_correctForIS) {
			recalculatepACoeffs(ionicStrength);
			ionicStrengthUnstable = !hasIonicStrengthConverged(X10(maxChargepACoeff), X10(m_pACoeffs.back()));
		} else
			ionicStrengthUnstable = false;

		ECHMET_DEBUG_CODE(std::cerr << m_pCx << "\n");

		ECHMET_DEBUG_CODE(fprintf(stderr, "Solver finished: Ionic strength = %g [mol/dm3]\n", CAESRealToDouble(ionicStrength)));


	} while (ionicStrengthUnstable);

	ECHMET_DEBUG_CODE(fprintf(stderr, "Total iterations = %u\n", totalIterationsCtr));

	m_finalIonicStrength = ionicStrength;

	m_totalIterations = totalIterationsCtr;
	m_outerIterations = outerIterationsCtr;

	return RetCode::OK;
}

/*!
 * Returns solved ionic concentrations as a matrix of \p CAESReal s
 *
 * @return A \p SolverVector<CAESReal> object.
 */
template <typename CAESReal, InstructionSet ISet>
SolverVector<CAESReal> SolverInternal<CAESReal, ISet>::rawConcentrations() const
{
	SolverVector<CAESReal> concentrations(m_pCx.rows());

	for (int rowCounter = 0; rowCounter < m_pCx.rows(); rowCounter++)
		concentrations(rowCounter) = X10(CVI(m_pCx, rowCounter));

	return concentrations;
}

/*!
 * Retunrs value of ionic strength as \p CAESReal.
 *
 * @return Value of ionic strength in <tt>mol/dm3</tt>.
 */
template <typename CAESReal, InstructionSet ISet>
CAESReal SolverInternal<CAESReal, ISet>::rawIonicStrength() const
{
	return m_finalIonicStrength;
}

/*!
 * Converts the raw results to \p SysComp::CalculatedProperties data.
 *
 * @param[in,out] calcProps The \p SysComp::CalculatedProperties struct to fill.
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::resultsToOutput(SysComp::CalculatedProperties &calcProps) noexcept
{
	/* This relies on the assumption that H+ and OH- are always on the top
	 * of the ionic concentrations vector */
	for (int rowCounter = 0; rowCounter < m_pCx.rows(); rowCounter++)
		(*calcProps.ionicConcentrations)[rowCounter] = CAESRealToECHMETReal(X10(CVI(m_pCx, rowCounter))); /* Concentrations of all constituents */

	calcProps.ionicStrength = CAESRealToECHMETReal(m_finalIonicStrength);
}

/*!
 * Checks if a matrix contains valid numbers
 *
 * Checks if a matrix contains valid numbers. Anything besides INF and NAN
 * is considered valid. An exception is thrown if an invalid value is found.
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::validateMatrix(const SolverMatrix<CAESReal> &m)
{
	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			if (VMath::isinf(m(row, col)))
				throw NumericErrorException(NumericErrorException::ExType::ET_INF);

			if (VMath::isnan(m(row, col)))
				throw NumericErrorException(NumericErrorException::ExType::ET_NAN);
		}
	}
}

/*!
 * Checks if a vector contains valid numbers
 *
 * Checks if a vector contains valid numbers. Anything besides INF and NAN
 * is considered valid. An exception is thrown if an invalid value is found.
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::validateVector(const typename NewtonRaphson<CAESReal, ISet>::TX &v)
{
	for (int row = 0; row < v.rows(); row++) {
		if (VMath::isinf(v(row)))
			throw NumericErrorException(NumericErrorException::ExType::ET_INF);

		if (VMath::isnan(v(row)))
			throw NumericErrorException(NumericErrorException::ExType::ET_NAN);
	}
}

/*!
 * Checks if a vector contains valid numbers
 *
 * Checks if a vector contains valid numbers. Anything besides INF and NAN
 * is considered valid. An exception is thrown if an invalid value is found.
 */
template <typename CAESReal, InstructionSet ISet>
void SolverInternal<CAESReal, ISet>::validateVector(const SolverVector<CAESReal> &v)
{
	for (int row = 0; row < v.rows(); row++) {
		if (VMath::isinf(v(row)))
			throw NumericErrorException(NumericErrorException::ExType::ET_INF);

		if (VMath::isnan(v(row)))
			throw NumericErrorException(NumericErrorException::ExType::ET_NAN);
	}
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_SOLVERINTERNAL_HPP
