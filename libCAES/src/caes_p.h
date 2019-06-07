#ifndef ECHMET_CAES_CAES_P_H
#define ECHMET_CAES_CAES_P_H

#include <echmetcaes.h>
#include "mappedmatrix.h"
#include "totalequilibrium.h"
#include "chargesummer.h"
#include "vecmath/vecmath.h"

#define ECHMET_IMPORT_INTERNAL
#include <mpreal.h>
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

#include <memory>

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class SolverContextImpl;

template <bool>
class FreeMPFRCacheSwitch {};

template <typename CAESReal, InstructionSet ISet>
class SolverInternal;

template <typename CAESReal, InstructionSet ISet>
inline
void releaseRawArray(CAESReal *v) noexcept
{
	AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(v);
}

template <typename CAESReal, InstructionSet ISet>
using RawArrayPtr = std::unique_ptr<CAESReal, decltype(&releaseRawArray<CAESReal, ISet>)>;

template <typename CAESReal, InstructionSet ISet>
RawArrayPtr<CAESReal, ISet> makeRawArray(size_t size)
{
	return RawArrayPtr<CAESReal, ISet>{AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(size), &releaseRawArray<CAESReal, ISet>};
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class UnsafeContext;

template <typename CAESReal, InstructionSet ISet>
class UnsafeContext<CAESReal, ISet, false> {
public:
	UnsafeContext();

	SolverInternal<CAESReal, ISet> *internal;
	SolverVector<CAESReal> *anCVec;

	RawArrayPtr<CAESReal, ISet> estimatedConcentrations;
	RawArrayPtr<CAESReal, ISet> estimatedIC;
	RawArrayPtr<CAESReal, ISet> dEstimatedICdH;

	std::vector<CAESReal> activityCoefficients;
	ChargeSummer<CAESReal, ISet, false> *chargeSummer;
};

template <typename CAESReal, InstructionSet ISet>
class UnsafeContext<CAESReal, ISet, true> {
public:
	/* Empty */
};

/*!
 * Equilibria solver object implementation
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class SolverImpl : public Solver {
public:
	explicit SolverImpl(SolverContextImpl<CAESReal> *ctx, const Options options, const NonidealityCorrections corrections);
	virtual ~SolverImpl() noexcept override;
	virtual SolverContext * ECHMET_CC context() noexcept override;
	SolverContextImpl<CAESReal> * contextInternal() noexcept;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual RetCode ECHMET_CC estimateDistributionFast(const ECHMETReal &cHInitial, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept override;
	virtual RetCode ECHMET_CC estimateDistributionSafe(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept override;
	template <typename OutputReal>
	RetCode fillResults(const std::pair<CAESReal *, CAESReal> &results,
			    OutputReal *estimatedConcentrations, OutputReal &ionicStrength) noexcept;
	virtual Options ECHMET_CC options() const noexcept override;
	virtual RetCode ECHMET_CC setContext(SolverContext *ctx) noexcept override;
	virtual RetCode ECHMET_CC setOptions(const Options options) noexcept override;
	virtual RetCode ECHMET_CC solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept override;
	RetCode solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const Vec<CAESReal> *estimatedConcentrations,
			 const size_t iterations, SolverIterations *iterationsNeeded = nullptr) noexcept;

	template <typename OutputReal>
	RetCode estimateDistributionSafeInternal(const RealVec *const ECHMET_RESTRICT_PTR analyticalConcentrations,
						 OutputReal *const ECHMET_RESTRICT_PTR estimatedConcentrations, OutputReal &ionicStrength) noexcept;

private:
	void defaultActivityCoefficients(std::vector<CAESReal> &activityCoefficients) const;

	std::pair<CAESReal *, CAESReal> estimatepHFast(const CAESReal &cHInitial, const ECHMETReal *analyticalConcentrations,
						       CAESReal *const ECHMET_RESTRICT_PTR icConcs, CAESReal *const ECHMET_RESTRICT_PTR dIcConcsdH,
						       std::vector<CAESReal> &activityCoefficients,
						       ChargeSummer<CAESReal, ISet, ThreadSafe> &chargeSummer);
	std::pair<CAESReal *, CAESReal> estimatepHSafe(const ECHMETReal *analyticalConcentrations,
						       CAESReal *const ECHMET_RESTRICT_PTR icConcs,
						       std::vector<CAESReal> &activityCoefficients,
						       ChargeSummer<CAESReal, ISet, ThreadSafe> &chargeSummer);
	void initializeTotalEquilibria(const SolverContextImpl<CAESReal> *ctx);
	SolverInternal<CAESReal, ISet> * makeSolverInternal(const SolverContextImpl<CAESReal> *ctx) const;
	void releaseSolverInternal(SolverInternal<CAESReal, ISet> *internal, FreeMPFRCacheSwitch<true>) noexcept;
	void releaseSolverInternal(SolverInternal<CAESReal, ISet> *internal, FreeMPFRCacheSwitch<false>) noexcept;
	RetCode setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept;

	SolverContextImpl<CAESReal> *m_ctx;			/*!< Associated solver context */

	Solver::Options m_options;

	const bool m_correctDebyeHuckel;			/*!< Correct with Debye-Hückel */

	size_t m_TECount;					/*!< Number of ionic concentrations that can be estimated from G-polynomial */
	std::vector<TotalEquilibrium<CAESReal, ThreadSafe>> m_totalEquilibria;

	size_t m_totalLigandCopySize;

	UnsafeContext<CAESReal, ISet, ThreadSafe> m_unsafe;

	template <typename X, InstructionSet, bool> friend class SolverImplSpec;
};

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class SolverImplSpec;

template <typename CAESReal, InstructionSet ISet>
class SolverImplSpec<CAESReal, ISet, false>
{
public:
	using SI = SolverImpl<CAESReal, ISet, false>;

	SolverImplSpec() = delete;

	static std::pair<CAESReal *, CAESReal> estimatepHFastWrapper(SI *solver, const CAESReal &cHInitial, const RealVec *analyticalConcentrations);
	static std::pair<CAESReal *, CAESReal> estimatepHSafeWrapper(SI *solver, const RealVec *analyticalConcentrations);
	static void initializeEstimators(SI *solver);
	static void initializeUnsafe(SI *solver, const SolverContextImpl<CAESReal> *ctx);
	static void releaseUnsafe(SI *solver) noexcept;
	static void setContainersForSolve(SI *solver,
			SolverInternal<CAESReal, ISet> *&internal,
				   SolverVector<CAESReal> *&anCVec,
				   CAESReal *&estimatedConcentrations);
	static void setContainersForSolveRaw(SI *solver,
					     SolverInternal<CAESReal, ISet> *&internal,
					     CAESReal *&estimatedConcentrationsInternal);
};

template <typename CAESReal, InstructionSet ISet>
class SolverImplSpec<CAESReal, ISet, true>
{
public:
	using SI = SolverImpl<CAESReal, ISet, true>;

	SolverImplSpec() = delete;

	static std::pair<CAESReal *, CAESReal> estimatepHFastWrapper(SI *solver, const CAESReal &cHInitial, const RealVec *analyticalConcentrations);
	static std::pair<CAESReal *, CAESReal> estimatepHSafeWrapper(SI *solver, const RealVec *analyticalConcentrations);
	static void initializeEstimators(SI *solver);
	static void initializeUnsafe(SI *solver, const SolverContextImpl<CAESReal> *ctx);
	static void releaseUnsafe(SI *solver) noexcept;
	static void setContainersForSolve(SI *solver,
			SolverInternal<CAESReal, ISet> *&internal,
				   SolverVector<CAESReal> *&anCVec,
				   CAESReal *&estimatedConcentrations);
	static void setContainersForSolveRaw(SI *solver,
					     SolverInternal<CAESReal, ISet> *&internal,
					     CAESReal *&estimatedConcentrationsInternal);
};

template <typename CAESReal, bool ThreadSafe>
static
void calculateDistribution(const CAESReal &v, SolverVector<CAESReal> &distribution, std::vector<TotalEquilibrium<CAESReal, ThreadSafe>> &totalEquilibria, const RealVec *analyticalConcentrations);

template <typename CAESReal>
static
void calculateMaximumVariants(const size_t i, const LigandIonicFormVec<CAESReal> &ligandIFs, uint32_t &total, uint32_t accum, const bool exclusive) noexcept;

template <typename CAESReal>
static
Solver * createSolverInternal(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) noexcept;

template <typename CAESReal>
static
RetCode createSolverContextInternal(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept;

InstructionSet detectInstructionSet() noexcept;

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
