#ifndef ECHMET_CAES_CAES_P_H
#define ECHMET_CAES_CAES_P_H

#include <echmetcaes.h>
#include <vector>
#include "types.h"

#define ECHMET_IMPORT_INTERNAL
#include <mpreal.h>
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class SolverContextImpl;

template <typename CAESReal>
class SolverInternal;

/*!
 * Description of total distribution equilibria for a given constituent
 */
template <typename CAESReal>
class TotalEquilibrium {
public:
	TotalEquilibrium();

	/*!
	 * TotalEquilibrium c-tor
	 *
	 * @param[in] numLow Lowest equilibrium index
	 * @param[in] numHigh Highest equilibirium index
	 * @param[in] pBs Vector of consecutive equilibrium constants
	 * @param[in] concentration Analytical concentration of the constituent
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<ECHMETReal> &pBs, const CAESReal &concentration) :
		concentration(concentration),
		Ls(calculateLs(pBs)),
		numLow(numLow),
		numHigh(numHigh)
	{
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const CAESReal &concentration = 1.0);
	TotalEquilibrium(TotalEquilibrium &&other);
	std::vector<CAESReal> concentrations(const CAESReal &v) const;
	std::vector<CAESReal> distribution(const CAESReal &v) const;
	std::vector<CAESReal> Ts(const CAESReal &v, CAESReal &X) const;

	const CAESReal concentration;		/*!< Analytical concentration of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */

private:
	/*!
	 * Calculates total equilibrium constants from constecutive constants
	 *
	 * @param[in] pBs Consecutive equilibrium constants
	 *
	 * @return Vector of total equilibrium constants
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	std::vector<CAESReal> calculateLs(const std::vector<ECHMETReal> &pBs)
	{
		std::vector<CAESReal> _pBs;
		_pBs.reserve(pBs.size());

		for (const ECHMETReal &d : pBs)
			_pBs.emplace_back(d);

		return calculateLs<CAESReal>(_pBs);
	}

	std::vector<CAESReal> calculateLs(const std::vector<CAESReal> &pBs);
};


template <bool>
class FreeMPFRCacheSwitch {};

/*!
 * Equilibria solver object implementation
 */
template <typename CAESReal>
class SolverImpl : public Solver {
public:
	explicit SolverImpl(const SolverContextImpl<CAESReal> *ctx, const NonidealityCorrections corrections, const Options options) noexcept;
	virtual ~SolverImpl() noexcept override;
	virtual const SolverContext * ECHMET_CC context() const noexcept override;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual Options ECHMET_CC options() const noexcept override;
	virtual RetCode ECHMET_CC setContext(const SolverContext *ctx) noexcept override;
	virtual RetCode ECHMET_CC setOptions(const Options options) noexcept override;
	virtual RetCode ECHMET_CC solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept override;
	RetCode solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const SolverVector<CAESReal> &estimatedConcentrations, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept;

private:
	void releaseSolverInternal(SolverInternal<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept;
	void releaseSolverInternal(SolverInternal<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept;

	Options m_options;					/*!< Solver options */
	bool m_correctDebyeHuckel;				/*!< Correct with Debye-Hückel */
	const SolverContextImpl<CAESReal> *m_ctx;		/*!< Associated solver context */
	SolverInternal<CAESReal> *m_internalUnsafe;		/*!< Internal solver used by thread-unsafe variant of the solver */

	SolverVector<CAESReal> m_anCVec;
	SolverVector<CAESReal> m_estimatedConcentrations;

	RetCode setContextInternal(const SolverContextImpl<CAESReal> *ctx) noexcept;
};

template <typename CAESReal>
void calculateDistribution(const CAESReal &v, SolverVector<CAESReal> &distribution, const std::vector<TotalEquilibrium<CAESReal>> &totalEquilibria);

template <typename CAESReal>
void calculateMaximumVariants(const size_t i, const LigandIonicFormVec<CAESReal> &ligandIFs, uint32_t &total, uint32_t accum, const bool exclusive) noexcept;

template <typename CAESReal>
Solver * createSolverInternal(const SolverContext *ctx, const Solver::Options options) noexcept;

template <typename CAESReal>
RetCode createSolverContextInternal(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept;

template <typename CAESReal>
void estimateComplexesDistribution(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands,
				   const SolverVector<CAESReal> &estConcentrations, const size_t LGBlockOffset, SysComp::IonicConcentrationVec *ionicConcentrations);
template <typename CAESReal>
RetCode estimateDistributionInternal(const SolverContext *ctx, RealVec *analyticalConcentrations, SolverVector<CAESReal> &estimatedConcentrations) noexcept;

template <typename CAESReal>
SolverVector<CAESReal> estimatepH(const std::vector<TotalEquilibrium<CAESReal>> &totalEquilibria);

template <typename CAESReal>
void generateComplexForms(Form<CAESReal> *f, FormVec<CAESReal> &forms, std::vector<FormVec<CAESReal>> &blocks, const LigandIonicFormVec<CAESReal> &ligandIFs, size_t gidx, const size_t stop);

template <typename CAESReal>
RetCode globalDataToInternal(LigandVec<CAESReal> *allLigands, LigandIonicFormVec<CAESReal> *allLigandIFs, CNVec<CAESReal> *cnVec, FormVec<CAESReal> *allForms,
			     const SysComp::IonicFormVec *gIfVec, const SysComp::ConstituentVec *gcVec) noexcept;
template <typename CAESReal>
RetCode initializeForms(ComplexNucleus<CAESReal> *cn, FormVec<CAESReal> *allForms, const SysComp::Constituent *ic, const SysComp::IonicFormVec *gIfVec, const LigandIonicFormVec<CAESReal> *allLigandIFs);

template <typename CAESReal>
SolverMatrix<CAESReal> * prepareJacobian(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands,
					 const size_t allFormsCount, const size_t allLigandIFsCount);

} // namespace CAES
} // namespace ECHMET

#include "caes.hpp"
#include "solverimpl.hpp"

#endif // ECHMET_CAES_CAES_P_H
