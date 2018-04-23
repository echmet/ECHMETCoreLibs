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
	explicit SolverImpl(SolverContextImpl<CAESReal> *ctx, const NonidealityCorrections corrections);
	virtual ~SolverImpl() noexcept override;
	virtual SolverContext * ECHMET_CC context() noexcept override;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual RetCode ECHMET_CC setContext(SolverContext *ctx) noexcept override;
	virtual RetCode ECHMET_CC solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept override;
	RetCode solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const SolverVector<CAESReal> &estimatedConcentrations, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept;

private:
	SolverInternalBase<CAESReal> * makeSolverInternal(const SolverContextImpl<CAESReal> *ctx) const;
	void releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept;
	void releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept;

	SolverContextImpl<CAESReal> *m_ctx;			/*!< Associated solver context */

	SolverInternalBase<CAESReal> *m_internalUnsafe;		/*!< Internal solver used by thread-unsafe variant of the solver */
	SolverVector<CAESReal> *m_anCVecUnsafe;			/*!< Vector of analytical concentrations used by thread-unsafe variant of the solver */
	CAESReal *m_estimatedConcentrationsUnsafe;		/*!< Aligned array of estimated concentrations used by thread-unsafe variant of the solver */

	const InstructionSet m_instructionSet;			/*!< Highest available CPU SIMD instruction set */
	bool m_correctDebyeHuckel;				/*!< Correct with Debye-Hückel */

	RetCode setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept;
};

template <typename CAESReal, bool ThreadSafe>
void calculateDistribution(const CAESReal &v, SolverVector<CAESReal> &distribution, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations);

template <typename CAESReal>
void calculateMaximumVariants(const size_t i, const LigandIonicFormVec<CAESReal> &ligandIFs, uint32_t &total, uint32_t accum, const bool exclusive) noexcept;

template <typename CAESReal>
Solver * createSolverInternal(SolverContext *ctx) noexcept;

template <typename CAESReal>
RetCode createSolverContextInternal(SolverContext *&ctx, const SolverContext::Options options, const SysComp::ChemicalSystem &chemSystem) noexcept;

template <typename CAESReal>
void estimateComplexesDistribution(const CNVec<CAESReal> *complexNuclei, const LigandVec<CAESReal> *allLigands,
				   const SolverVector<CAESReal> &estConcentrations, const size_t LGBlockOffset, SysComp::IonicConcentrationVec *ionicConcentrations);
template <typename CAESReal>
RetCode estimateDistributionInternal(const CAESReal &cHInitial, SolverContext *ctx, RealVec *analyticalConcentrations, SolverVector<CAESReal> &estimatedConcentrations,
				     const bool useFastEstimate) noexcept;

template <typename CAESReal, bool ThreadSafe>
SolverVector<CAESReal> estimatepHFast(const CAESReal &cHInitial, std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations,
				      SolverVector<CAESReal> &icConcs, SolverVector<CAESReal> &dIcConcsdH);

template <typename CAESReal, bool ThreadSafe>
SolverVector<CAESReal> estimatepHSafe(std::vector<TotalEquilibriumBase *> &totalEquilibria, const RealVec *analyticalConcentrations, SolverVector<CAESReal> &icConcs);

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
