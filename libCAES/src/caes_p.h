#ifndef ECHMET_CAES_CAES_P_H
#define ECHMET_CAES_CAES_P_H

#include <echmetcaes.h>
#include "types.h"
#include "totalequilibrium.h"
#include "vecmath/vecmath.h"

#define ECHMET_IMPORT_INTERNAL
#include <mpreal.h>
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class SolverContextImpl;

template <bool>
class FreeMPFRCacheSwitch {};

template <typename CAESReal>
class SolverInternalBase;

/*!
 * Equilibria solver object implementation
 */
template <typename CAESReal>
class SolverImpl : public Solver {
public:
	explicit SolverImpl(SolverContextImpl<CAESReal> *ctx, const Options options, const NonidealityCorrections corrections);
	virtual ~SolverImpl() noexcept override;
	virtual SolverContext * ECHMET_CC context() noexcept override;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual RetCode ECHMET_CC estimateDistributionFast(const ECHMETReal &cHInitial, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept override;
	virtual RetCode ECHMET_CC estimateDistributionSafe(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept override;
	RetCode estimateDistributionInternal(const CAESReal &cHInitial, const RealVec *analyticalConcentrations, SolverVector<CAESReal> &estimatedConcentrations,
					     CAESReal &ionicStrength, const bool useFastEstimate) noexcept;
	virtual Options ECHMET_CC options() const noexcept override;
	virtual RetCode ECHMET_CC setContext(SolverContext *ctx) noexcept override;
	virtual RetCode ECHMET_CC setOptions(const Options options) noexcept override;
	virtual RetCode ECHMET_CC solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept override;
	RetCode solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const SolverVector<CAESReal> &estimatedConcentrations,
			 const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept;

private:
	void defaultActivityCoefficients(std::vector<CAESReal> &activityCoefficients) const;
	template <bool ThreadSafe>
	std::pair<SolverVector<CAESReal>, CAESReal> estimatepHFast(const CAESReal &cHInitial, const RealVec *analyticalConcentrations,
								   SolverVector<CAESReal> &icConcs, SolverVector<CAESReal> &dIcConcsdH,
								   std::vector<CAESReal> &activityCoefficients);
	template <bool ThreadSafe>
	std::pair<SolverVector<CAESReal>, CAESReal> estimatepHSafe(const RealVec *analyticalConcentrations, SolverVector<CAESReal> &icConcs, std::vector<CAESReal> &activityCoefficients);
	void initializeEstimators();
	SolverInternalBase<CAESReal> * makeSolverInternal(const SolverContextImpl<CAESReal> *ctx) const;
	void initializeTotalEquilibria(const SolverContextImpl<CAESReal> *ctx);
	void releaseTotalEquilibria();
	void releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept;
	void releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept;
	RetCode setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept;

	SolverContextImpl<CAESReal> *m_ctx;			/*!< Associated solver context */

	Solver::Options m_options;

	SolverInternalBase<CAESReal> *m_internalUnsafe;		/*!< Internal solver used by thread-unsafe variant of the solver */
	SolverVector<CAESReal> *m_anCVecUnsafe;			/*!< Vector of analytical concentrations used by thread-unsafe variant of the solver */
	CAESReal *m_estimatedConcentrationsUnsafe;		/*!< Aligned array of estimated concentrations used by thread-unsafe variant of the solver */
	const bool m_correctDebyeHuckel;			/*!< Correct with Debye-Hückel */

	const InstructionSet m_instructionSet;			/*!< Highest available CPU SIMD instruction set */

	size_t m_TECount;					/*!< Number of ionic concentrations that can be estimated from G-polynomial */
	std::vector<TotalEquilibriumBase *> m_totalEquilibria;
	SolverVector<CAESReal> m_estimatedIonicConcentrations;
	SolverVector<CAESReal> m_dEstimatedIonicConcentrationsdH;
	std::vector<CAESReal> m_activityCoefficients;
};

template <typename CAESReal, bool ThreadSafe>
static
void calculateDistribution(const CAESReal &v, SolverVector<CAESReal> &distribution, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations);

template <typename CAESReal>
static
void calculateMaximumVariants(const size_t i, const LigandIonicFormVec<CAESReal> &ligandIFs, uint32_t &total, uint32_t accum, const bool exclusive) noexcept;

template <typename CAESReal>
static
Solver * createSolverInternal(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) noexcept;

template <typename CAESReal>
static
RetCode createSolverContextInternal(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept;

template <typename CAESReal>
static
void generateComplexForms(Form<CAESReal> *f, FormVec<CAESReal> &forms, std::vector<FormVec<CAESReal>> &blocks, const LigandIonicFormVec<CAESReal> &ligandIFs, size_t gidx, const size_t stop);

template <typename CAESReal>
static
RetCode globalDataToInternal(LigandVec<CAESReal> *allLigands, LigandIonicFormVec<CAESReal> *allLigandIFs, CNVec<CAESReal> *cnVec, FormVec<CAESReal> *allForms,
			     const SysComp::IonicFormVec *gIfVec, const SysComp::ConstituentVec *gcVec) noexcept;
template <typename CAESReal>
static
RetCode initializeForms(ComplexNucleus<CAESReal> *cn, FormVec<CAESReal> *allForms, const SysComp::Constituent *ic, const SysComp::IonicFormVec *gIfVec, const LigandIonicFormVec<CAESReal> *allLigandIFs);

template <typename CAESReal>
static
SolverMatrix<CAESReal> * prepareJacobian(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands,
					 const size_t allFormsCount, const size_t allLigandIFsCount);

} // namespace CAES
} // namespace ECHMET

#include "caes.hpp"
#include "solverimpl.hpp"

#endif // ECHMET_CAES_CAES_P_H
