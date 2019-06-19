#ifndef ECHMET_CAES_SOLVERIMPL_HPP
#define ECHMET_CAES_SOLVERIMPL_HPP

#include "solverinternal.h"
#include "estimator_helpers.hpp"
#include <functional>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	#include <x86intrin.h>
#else
	#include <xmmintrin.h>
	#include <immintrin.h>
#endif // ECHMET_COMPILER_

namespace ECHMET {
namespace CAES {

class FastEstimateFailureException : public std::exception {
public:
	using std::exception::exception;
};

template <typename CAESReal, bool V = std::is_same<CAESReal, mpfr::mpreal>::value>
static
FreeMPFRCacheSwitch<V> freeMPFRCache()
{
	return FreeMPFRCacheSwitch<V>{};
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void SolverImpl<CAESReal, ISet, ThreadSafe>::releaseSolverInternal(SolverInternal<CAESReal, ISet> *internal, FreeMPFRCacheSwitch<true>) noexcept
{
	if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY)) {
		delete internal;
		mpfr_free_cache();
	}
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void SolverImpl<CAESReal, ISet, ThreadSafe>::releaseSolverInternal(SolverInternal<CAESReal, ISet> *internal, FreeMPFRCacheSwitch<false>) noexcept
{
	if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY))
		delete internal;
}

template <typename CAESReal, InstructionSet ISet>
UnsafeContext<CAESReal, ISet, false>::UnsafeContext() :
	internal{nullptr},
	anCVec{nullptr},
	estimatedConcentrations{nullptr, releaseRawArray<CAESReal, ISet>},
	estimatedIC{nullptr, releaseRawArray<CAESReal, ISet>},
	dEstimatedICdH{nullptr, releaseRawArray<CAESReal, ISet>},
	chargeSummer{nullptr}
{
}

template <typename CAESReal, InstructionSet ISet>
std::pair<CAESReal *, CAESReal>
SolverImplSpec<CAESReal, ISet, false>::estimatepHFastWrapper(SolverImplSpec::SI *solver, const CAESReal &cHInitial, const RealVec *analyticalConcentrations)
{
	return solver->estimatepHFast(cHInitial, analyticalConcentrations->cdata(),
				      solver->m_unsafe.estimatedIC.get(), solver->m_unsafe.dEstimatedICdH.get(),
				      solver->m_unsafe.activityCoefficients, *solver->m_unsafe.chargeSummer);
}

template <typename CAESReal, InstructionSet ISet>
std::pair<CAESReal *, CAESReal>
SolverImplSpec<CAESReal, ISet, true>::estimatepHFastWrapper(SolverImplSpec::SI *solver, const CAESReal &cHInitial, const RealVec *analyticalConcentrations)
{
	std::vector<CAESReal> activityCoefficients;
	auto icConcs = makeRawArray<CAESReal, ISet>(solver->m_TECount);
	auto dIcConcsdH = makeRawArray<CAESReal, ISet>(solver->m_TECount);
	ChargeSummer<CAESReal, ISet, true> chargeSummer{solver->m_TECount, solver->m_totalEquilibria};

	activityCoefficients.resize(solver->m_ctx->chargesSquared.size());

	const auto ret = solver->estimatepHFast(cHInitial, analyticalConcentrations->cdata(),
						icConcs.get(), dIcConcsdH.get(),
						activityCoefficients, chargeSummer);
	icConcs.release();

	return ret;
}

template <typename CAESReal, InstructionSet ISet>
std::pair<CAESReal *, CAESReal>
SolverImplSpec<CAESReal, ISet, false>::estimatepHSafeWrapper(SolverImplSpec::SI *solver, const RealVec *analyticalConcentrations)
{
	return solver->estimatepHSafe(analyticalConcentrations->cdata(), solver->m_unsafe.estimatedIC.get(),
				      solver->m_unsafe.activityCoefficients, *solver->m_unsafe.chargeSummer);
}

template <typename CAESReal, InstructionSet ISet>
std::pair<CAESReal *, CAESReal>
SolverImplSpec<CAESReal, ISet, true>::estimatepHSafeWrapper(SolverImplSpec::SI *solver, const RealVec *analyticalConcentrations)
{
	std::vector<CAESReal> activityCoefficients;
	auto icConcs = makeRawArray<CAESReal, ISet>(solver->m_TECount);
	ChargeSummer<CAESReal, ISet, true> chargeSummer{solver->m_TECount, solver->m_totalEquilibria};

	activityCoefficients.resize(solver->m_ctx->chargesSquared.size());

	const auto ret = solver->estimatepHSafe(analyticalConcentrations->cdata(),
						icConcs.get(), activityCoefficients, chargeSummer);
	icConcs.release();

	return ret;
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, false>::initializeEstimators(SolverImplSpec::SI *solver)
{
	solver->m_unsafe.estimatedIC = makeRawArray<CAESReal, ISet>(solver->m_TECount);
	solver->m_unsafe.dEstimatedICdH = makeRawArray<CAESReal, ISet>(solver->m_TECount);

	solver->m_unsafe.chargeSummer = new ChargeSummer<CAESReal, ISet, false>{solver->m_TECount, solver->m_totalEquilibria};

	solver->m_unsafe.activityCoefficients.resize(solver->m_ctx->chargesSquared.size());
	solver->defaultActivityCoefficients(solver->m_unsafe.activityCoefficients);

	solver->m_totalLigandCopySize = 0;
	for (const auto l : *solver->m_ctx->allLigands)
		solver->m_totalLigandCopySize += l->chargeHigh - l->chargeLow + 1;
	solver->m_totalLigandCopySize *= sizeof(double);
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, true>::initializeEstimators(SolverImplSpec::SI *solver)
{
	solver->m_totalLigandCopySize = 0;
	for (const auto l : *solver->m_ctx->allLigands)
		solver->m_totalLigandCopySize += l->chargeHigh - l->chargeLow + 1;
	solver->m_totalLigandCopySize *= sizeof(double);
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, false>::initializeUnsafe(SolverImplSpec::SI *solver, const SolverContextImpl<CAESReal> *ctx)
{
	initializeEstimators(solver);

	solver->m_unsafe.internal = solver->makeSolverInternal(ctx);
	solver->m_unsafe.anCVec = new SolverVector<CAESReal>(ctx->analyticalConcentrationCount);
	solver->m_unsafe.estimatedConcentrations = makeRawArray<CAESReal, ISet>(ctx->concentrationCount);
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, true>::initializeUnsafe(SolverImplSpec::SI *, const SolverContextImpl<CAESReal> *)
{
	/* NOOP */
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, false>::releaseUnsafe(SolverImplSpec::SI *solver) noexcept
{
	delete solver->m_unsafe.internal;
	delete solver->m_unsafe.anCVec;

	solver->m_unsafe.internal = nullptr;
	solver->m_unsafe.anCVec = nullptr;
	solver->m_unsafe.estimatedConcentrations = nullptr;
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, true>::releaseUnsafe(SolverImplSpec::SI *) noexcept
{
	/* NOOP */
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, false>::setContainersForSolve(SolverImplSpec::SI *solver,
		SolverInternal<CAESReal, ISet> *&internal,
							      SolverVector<CAESReal> *&anCVec,
							      CAESReal *&estimatedConcentrations)
{
	internal = solver->m_unsafe.internal;
	anCVec = solver->m_unsafe.anCVec;
	estimatedConcentrations = solver->m_unsafe.estimatedConcentrations.get();
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, true>::setContainersForSolve(SolverImplSpec::SI *solver,
		SolverInternal<CAESReal, ISet> *&internal,
							     SolverVector<CAESReal> *&anCVec,
							     CAESReal *&estimatedConcentrations)
{
	internal = solver->makeSolverInternal(solver->m_ctx);
	anCVec = new SolverVector<CAESReal>(solver->m_ctx->analyticalConcentrationCount);
	estimatedConcentrations = AlignedAllocator<CAESReal, 32>::alloc(solver->m_ctx->concentrationCount);
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, false>::setContainersForSolveRaw(SI *solver,
								     SolverInternal<CAESReal, ISet> *&internal,
								     CAESReal *&estimatedConcentrationsInternal)
{
	internal = solver->m_unsafe.internal;
	estimatedConcentrationsInternal = solver->m_unsafe.estimatedConcentrations;
}

template <typename CAESReal, InstructionSet ISet>
void SolverImplSpec<CAESReal, ISet, true>::setContainersForSolveRaw(SI *solver,
								    SolverInternal<CAESReal, ISet> *&internal,
								    CAESReal *&estimatedConcentrationsInternal)
{
	internal = solver->makeSolverInternal(solver->m_ctx);
	estimatedConcentrationsInternal = AlignedAllocator<CAESReal, 32>::alloc(solver->m_ctx->concentrationCount);
}

/*!
 * Solver c-tor.
 *
 * @param[in] ctx \p SolverContext used to initialize solver's internal data structures. The pointer
 *            shall remain valid throughout the entire lifetime of the \p Solver .
 * @param[in] options Solver options.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
SolverImpl<CAESReal, ISet, ThreadSafe>::SolverImpl(SolverContextImpl<CAESReal> *ctx, const Options options, const NonidealityCorrections corrections) :
	m_ctx(ctx),
	m_options(options),
	m_correctDebyeHuckel(nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL))
{
	initializeTotalEquilibria(ctx);

	try {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::initializeUnsafe(this, ctx);
	} catch (const std::bad_alloc &) {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::releaseUnsafe(this);

		throw;
	}
}

/*!
 * Solver d-tor.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
SolverImpl<CAESReal, ISet, ThreadSafe>::~SolverImpl() ECHMET_NOEXCEPT
{
	SolverImplSpec<CAESReal, ISet, ThreadSafe>::releaseUnsafe(this);
}

/*!
 * Returns a pointer to \p SolverContext assigned to the solve.
 *
 * @return Pointer to \p SolverContext object.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
SolverContext * ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::context() ECHMET_NOEXCEPT
{
	return m_ctx;
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
SolverContextImpl<CAESReal> * SolverImpl<CAESReal, ISet, ThreadSafe>::contextInternal() ECHMET_NOEXCEPT
{
	return m_ctx;
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::destroy() const ECHMET_NOEXCEPT
{
	delete this;
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void SolverImpl<CAESReal, ISet, ThreadSafe>::defaultActivityCoefficients(std::vector<CAESReal> &activityCoefficients) const
{
	std::fill(activityCoefficients.begin(), activityCoefficients.end(), 1.0);
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::estimateDistributionFast(const ECHMETReal &cHInitial, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	if (analyticalConcentrations->size() != m_ctx->analyticalConcentrationCount)
		return RetCode::E_INVALID_ARGUMENT;

	assert(calcProps.ionicConcentrations->size() == m_ctx->concentrationCount);

	try {
		const auto results = SolverImplSpec<CAESReal, ISet, ThreadSafe>::estimatepHFastWrapper(this, cHInitial, analyticalConcentrations);

		fillResults<ECHMETReal>(results, calcProps.ionicConcentrations->data(), calcProps.ionicStrength);
	} catch (const FastEstimateFailureException &) {
		return RetCode::E_FAST_ESTIMATE_FAILURE;
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	return RetCode::OK;
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::estimateDistributionSafe(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	if (analyticalConcentrations->size() != m_ctx->analyticalConcentrationCount)
		return RetCode::E_INVALID_ARGUMENT;

	assert(calcProps.ionicConcentrations->size() == m_ctx->concentrationCount);

	return estimateDistributionSafeInternal<ECHMETReal>(analyticalConcentrations, calcProps.ionicConcentrations->data(), calcProps.ionicStrength);
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe> template <typename OutputReal>
RetCode SolverImpl<CAESReal, ISet, ThreadSafe>::estimateDistributionSafeInternal(const RealVec *const ECHMET_RESTRICT_PTR analyticalConcentrations,
										 OutputReal *const ECHMET_RESTRICT_PTR estimatedConcentrations,
										 OutputReal &ionicStrength) noexcept
{
	try {
		const auto results = SolverImplSpec<CAESReal, ISet, ThreadSafe>::estimatepHSafeWrapper(this, analyticalConcentrations);

		fillResults<OutputReal>(results, estimatedConcentrations, ionicStrength);
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	return RetCode::OK;
}

/*!
 * Calculates the initial estimation of concentration of all species in the system.
 *
 * @param[in] cHInitial Pre-estimated value of concentration of H<sub>3</sub>O<sup>+</sup> ions.
 *                      This is used only in fast estimate.
 * @param[in] analyticalConcentrations Vector of analytical concentrations of constituents.
 * @param[in, out] estimatedConcentrations Output vector of estimated ionic concentrations of all species.
 *                                         The vector will be resized appropriately by this function.
 * @param[in, out] ionicStrength Estimated value of ionic strength.
 * @param[in] useFastEstimate When <tt>true</tt>, fast estimate will be performed.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of the concentration vectors.
 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution.
 * @retval RetCode::E_FAST_ESTIMATE_FAILURE Fast estimation failed to find a solution.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe> template <typename OutputReal>
RetCode SolverImpl<CAESReal, ISet, ThreadSafe>::fillResults(const std::pair<CAESReal *, CAESReal> &results, OutputReal *estimatedConcentrations, OutputReal &ionicStrength) noexcept
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "Estimated pH = %g\n", CAESRealToDouble(pX(results.first[0]) + 3.0)));
	ECHMET_DEBUG_CODE(fprintf(stderr, "Estimated ionic strength = %g\n", CAESRealToDouble(ionicStrength)));

	/* H+ and OH- are expected to be the first and second item in the vector */
	estimatedConcentrations[0] = CAESRealToECHMETReal<CAESReal, OutputReal>(results.first[0]);
	estimatedConcentrations[1] = CAESRealToECHMETReal<CAESReal, OutputReal>(results.first[1]);

	estimateComplexesDistribution<CAESReal, OutputReal>(m_ctx->complexNuclei, m_ctx->allLigands, m_totalLigandCopySize,
							    results.first, m_ctx->allForms->size() + 2, estimatedConcentrations);

	ionicStrength = CAESRealToECHMETReal<CAESReal, OutputReal>(results.second);

	if (ThreadSafe)
		releaseRawArray<CAESReal, ISet>(results.first);

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
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
std::pair<CAESReal *, CAESReal> SolverImpl<CAESReal, ISet, ThreadSafe>::estimatepHFast(const CAESReal &cHInitial, const ECHMETReal *analyticalConcentrations,
										       CAESReal *const ECHMET_RESTRICT_PTR icConcs, CAESReal *const ECHMET_RESTRICT_PTR dIcConcsdH,
										       std::vector<CAESReal> &activityCoefficients,
										       ChargeSummer<CAESReal, ISet, ThreadSafe> &chargeSummer)
{
	const CAESReal KW_298 = CAESReal(PhChConsts::KW_298) * 1e6;
	const CAESReal threshold = electroneturalityPrecision<CAESReal>();
	CAESReal cH = cHInitial;
	CAESReal cHNew;

	CAESReal ionicStrength = 0.0;
	size_t isLoopCtr = 0;
	CAESReal maxChargeActivityCoeff;
	bool ionicStrengthUnstable;

	if (m_correctDebyeHuckel || ThreadSafe)
		defaultActivityCoefficients(activityCoefficients);

	do {
		const CAESReal activityOneSquared = activityCoefficients[1] * activityCoefficients[1];
		size_t ctr = 0;
		maxChargeActivityCoeff = activityCoefficients.back();
		CAESReal cOH;
		CAESReal z;
		CAESReal dZ;

		while (true) {
			calculateDistributionWithDerivative<CAESReal, ISet, ThreadSafe>(cH, icConcs, dIcConcsdH,
											m_totalEquilibria,
											analyticalConcentrations,
											activityCoefficients);

			cOH = KW_298 / (cH * activityOneSquared);

			icConcs[0] = cH;
			icConcs[1] = cOH;
			dIcConcsdH[0] = 1.0;
			dIcConcsdH[1] = cOH * cOH;

			chargeSummer.calcWithdZ(icConcs, dIcConcsdH, z, dZ);

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

		ionicStrength = chargeSummer.calculateIonicStrength(icConcs);
		if (m_correctDebyeHuckel) {
			calculateActivityCoefficients(ionicStrength, activityCoefficients, m_ctx->chargesSquared);
			ionicStrengthUnstable = !hasIonicStrengthConverged(maxChargeActivityCoeff, activityCoefficients.back());
		} else
			ionicStrengthUnstable = false;
	} while (m_correctDebyeHuckel && ionicStrengthUnstable && isLoopCtr++ < 100);

	return { icConcs, ionicStrength };
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
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
std::pair<CAESReal *, CAESReal> SolverImpl<CAESReal, ISet, ThreadSafe>::estimatepHSafe(const ECHMETReal *analyticalConcentrations,
										       CAESReal *const ECHMET_RESTRICT_PTR icConcs,
										       std::vector<CAESReal> &activityCoefficients,
										       ChargeSummer<CAESReal, ISet, ThreadSafe> &chargeSummer)
{
	const CAESReal KW_298 = CAESReal(PhChConsts::KW_298) * 1e6;
	const CAESReal threshold = electroneturalityPrecision<CAESReal>();

	CAESReal ionicStrength = 0.0;
	size_t isLoopCtr = 0;
	CAESReal maxChargeActivityCoeff;
	bool ionicStrengthUnstable;

	if (m_correctDebyeHuckel || ThreadSafe)
		defaultActivityCoefficients(activityCoefficients);

	do {
		const CAESReal activityOneSquared = activityCoefficients[1] * activityCoefficients[1];
		CAESReal cH = 1.0e-4;
		CAESReal leftWall = 0.0;
		CAESReal rightWall = 100000.0;
		size_t ctr = 0;
		maxChargeActivityCoeff = activityCoefficients.back();

		CAESReal cOH;

		while (true) {
			calculateDistribution<CAESReal, ISet, ThreadSafe>(cH, icConcs, m_totalEquilibria, analyticalConcentrations,
									  activityCoefficients);

			cOH = KW_298 / (cH * activityOneSquared);

			icConcs[0] = cH;
			icConcs[1] = cOH;

			const CAESReal z = chargeSummer.calc(icConcs);

			//fprintf(stderr, "cH %g, z %g, dZ %g, cHNew, %g\n", CAESRealToDouble(cH), CAESRealToDouble(z), CAESRealToDouble(dZ), CAESRealToDouble(cHNew));

			if (VMath::abs(z) < threshold)
				break;

			/* Maximum number of iterations exceeded, return what we have so far and
			 * hope for the best */
			if (ctr++ > 1000)
				break;

			if (z > 0)
				rightWall = cH;
			else
				leftWall = cH;

			cH = (rightWall - leftWall) / 2.0 + leftWall;
		}

		ionicStrength = chargeSummer.calculateIonicStrength(icConcs);
		if (m_correctDebyeHuckel) {
			calculateActivityCoefficients(ionicStrength, activityCoefficients, m_ctx->chargesSquared);
			ionicStrengthUnstable = !hasIonicStrengthConverged(maxChargeActivityCoeff, activityCoefficients.back());
		} else
			ionicStrengthUnstable = false;
	} while (m_correctDebyeHuckel && ionicStrengthUnstable && isLoopCtr++ < 100);

	return { icConcs, ionicStrength };
}

/*!
 * Initializes TotalEquilibria objects used to estimate ionic distribution
 *
 * @param[in] ctx SovlerContextImpl to use
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
void SolverImpl<CAESReal, ISet, ThreadSafe>::initializeTotalEquilibria(const SolverContextImpl<CAESReal> *ctx)
{
	auto countTEForms = [](const auto &TEVec) {
		size_t TECount = 0;
		for (const auto &te : TEVec) {
			TECount += te.numHigh - te.numLow + 1;
		}
		return TECount;
	};

	m_totalEquilibria.clear();

	m_totalEquilibria.reserve(ctx->complexNuclei->size() + ctx->allLigands->size());

	for (const ComplexNucleus<CAESReal> *cn : *ctx->complexNuclei)
		m_totalEquilibria.emplace_back(cn->chargeLow, cn->chargeHigh, cn->pKas, cn->analyticalConcentrationIndex);
	for (const Ligand<CAESReal> *l : *ctx->allLigands)
		m_totalEquilibria.emplace_back(l->chargeLow, l->chargeHigh, l->pKas, l->analyticalConcentrationIndex);

	m_TECount = countTEForms(m_totalEquilibria);

	m_TECount += 2;
}

/*!
 * Makes appropriate \p SolverInternal object
 * based on available SIMD instruction set
 *
 * @retval Pointer to \p SolverInternal object
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
SolverInternal<CAESReal, ISet> * SolverImpl<CAESReal, ISet, ThreadSafe>::makeSolverInternal(const SolverContextImpl<CAESReal> *ctx) const
{
	return new SolverInternal<CAESReal, ISet>(ctx);
}

/*!
 * Returns the current options of the solver.
 *
 * @return Solver options.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
typename SolverImpl<CAESReal, ISet, ThreadSafe>::Options ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::options() const noexcept
{
	return m_options;
}

/*!
 * Sets new \p SolverContext.
 *
 * @param[in] ctx Pointer to the new \p SolverContext.
 * @param[in] totalConcentrations Pointer to vector of new total concentrations of all system constituents.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT SolverContext not castable to SolverContextImpl.
 * @retval RetCode::E_INVALID_ARGUMENT Change of thread safetiness is not allowed.
 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate internal resources.
 *				If this error code is returned the sovler must be assumed to be in
 *				an inconsistent state.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::setContext(SolverContext *ctx) ECHMET_NOEXCEPT
{
	SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	return setContextInternal(ctxImpl);
}

/*!
 * Sets new \p SolverContext. Internal implementation.
 * Failure of this function may leave the solver in an inconsistent state.
 *
 * @param[in] ctx Pointer to the new \p SolverContextImpl.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate internal resources.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode SolverImpl<CAESReal, ISet, ThreadSafe>::setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept
{
	m_totalEquilibria.clear();

	try {
		initializeTotalEquilibria(ctx);
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	SolverImplSpec<CAESReal, ISet, ThreadSafe>::releaseUnsafe(this);
	try {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::initializeUnsafe(this, ctx);
	} catch (const std::bad_alloc &) {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::releaseUnsafe(this);

		return RetCode::E_NO_MEMORY;
	}

	m_ctx = ctx;

	return RetCode::OK;
}

/*!
 * Sets new solver options.
 *
 * @param[in] options New options.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Nonsensical option value was passed as the argument.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::setOptions(const Options options) noexcept
{
	/* Changing thread-safetiness is not allowed */
	if ((m_options & Options::DISABLE_THREAD_SAFETY) != (options & Options::DISABLE_THREAD_SAFETY))
		return RetCode::E_INVALID_ARGUMENT;

	m_options = options;

	return RetCode::OK;
}

/*!
 * Calculates the equilibrium ionic distribution of the system. Internal implementation.
 *
 * @param[in] analyticalConcentrations Analytical concentrations of constituents in the system.
 * @param[in,out] calcProps Struct where to store the results.
 * @param[in] Maximum number of iterations to try.
 * @param[in,out] iterationsNeeded If given it returns the number of iterations needed for the solver to converge. The passed object is not altered if the solver fails to converge.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Not enough memory to perform the calculation.
 * @retval RetCode::E_NRS_FAILURE Newton-Rapshon solver encountered an error during calculation.
 * @retval RetCode::E_NRS_NO_CONVERGENCE Newton-Raphson solver failed to converge within the given number of iterations
 * @retval RetCode::E_NRS_STUCK Greatest change of X-value calculated by the Newton-Raphson solver is below the precision threshold.
 * @retval RetCode::E_NRS_NO_SOLUTION System appears to have no solution.
 * @retval RetCode::E_IS_NO_CONVERGENCE Solver failed to find a solution within the given number of iterations.
 * @retval RetCode::E_INVALID_ARGUMENT Vector of analytical concentrations has invalid size or initial value of ionic strength is invalid.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode ECHMET_CC SolverImpl<CAESReal, ISet, ThreadSafe>::solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps,
								const size_t iterations, SolverIterations *iterationsNeeded) ECHMET_NOEXCEPT
{
	SolverInternal<CAESReal, ISet> *internal = nullptr;
	SolverVector<CAESReal> *anCVec = nullptr;
	CAESReal *estimatedConcentrations = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (ThreadSafe) {
			delete anCVec;
			AlignedAllocator<CAESReal, 32>::free(estimatedConcentrations);
		}
	};

	try {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::setContainersForSolve(this, internal, anCVec, estimatedConcentrations);
	} catch (const std::bad_alloc &) {
		delete anCVec;
		delete internal;
		return RetCode::E_NO_MEMORY;
	}

	if (analyticalConcentrations->size() != static_cast<size_t>(anCVec->rows())) {
		releaseResources();

		return RetCode::E_INVALID_ARGUMENT;
	}

	try {
		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			(*anCVec)(idx) = analyticalConcentrations->elem(idx);

		for (size_t idx = 0; idx < calcProps.ionicConcentrations->size(); idx++)
			estimatedConcentrations[idx] = calcProps.ionicConcentrations->elem(idx);

		const RetCode tRet = internal->solve(anCVec, estimatedConcentrations, m_correctDebyeHuckel, iterations, calcProps.ionicStrength);
		if (tRet != RetCode::OK) {
			releaseResources();

			return tRet;
		}

		internal->resultsToOutput(calcProps);

		if (iterationsNeeded != nullptr)
			*iterationsNeeded = internal->iterations();

		releaseResources();

		return RetCode::OK;
	} catch (std::bad_alloc &) {
		releaseResources();

		return RetCode::E_NO_MEMORY;
	}
}

/*!
 * Calculates the equilibrium ionic distribution of the system. Internal implementation with raw output.
 *
 * @param[in,out] concentrations Vector where the resulting concentrations will be stored. The concentrations are sorted in the
 *		  same order as in the <tt>CalculatedProperties::ionicConcentrations</tt> vector.
 * @param[in, out] ionicStrength Ionic strength of the solved system. If correction for Debye-Hückel is enabled,
 *                 the input value shall be set to ionic strength of the estimated system.
 * @param[in] anCVec Analytical concentrations of constituents in the system.
 * @param[in,out] estimatedConcentrations Vector of estimated concentrations.
 * @param[in] iterations Maximum number of iterations to try.
 * @param[in,out] iterationsNeeded If given it returns the number of iterations needed for the solver to converge. The passed object is not altered if the solver fails to converge.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Not enough memory to perform the calculation.
 * @retval RetCode::E_NRS_FAILURE Newton-Rapshon solver encountered an error during calculation.
 * @retval RetCode::E_NRS_NO_CONVERGENCE Newton-Raphson solver failed to converge within the given number of iterations
 * @retval RetCode::E_NRS_STUCK Greatest change of X-value calculated by the Newton-Raphson solver is below the precision threshold.
 * @retval RetCode::E_NRS_NO_SOLUTION System appears to have no solution.
 * @retval RetCode::E_IS_NO_CONVERGENCE Solver failed to find a solution within the given number of iterations.
 * @retval RetCode::E_INVALID_ARGUMENT Initial value of ionic strength is invalid.
 */
template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
RetCode SolverImpl<CAESReal, ISet, ThreadSafe>::solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec,
							 const Vec<CAESReal> *estimatedConcentrations, const size_t iterations, SolverIterations *iterationsNeeded) noexcept
{
	SolverInternal<CAESReal, ISet> *internal = nullptr;
	CAESReal *estimatedConcentrationsInternal = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY))
			AlignedAllocator<CAESReal, 32>::free(estimatedConcentrationsInternal);
	};

	try {
		SolverImplSpec<CAESReal, ISet, ThreadSafe>::setContainersForSolveRaw(this, internal, estimatedConcentrationsInternal);
	} catch (const std::bad_alloc &) {
		delete internal;
		return RetCode::E_NO_MEMORY;
	}

	for (size_t idx = 0; idx < estimatedConcentrations->size(); idx++)
		estimatedConcentrationsInternal[idx] = estimatedConcentrations->elem(idx);

	const RetCode tRet = internal->solve(anCVec, estimatedConcentrationsInternal, m_correctDebyeHuckel, iterations, ionicStrength);
	if (tRet != RetCode::OK) {
		releaseResources();

		return tRet;
	}

	try {
		concentrations = internal->rawConcentrations();
	} catch (std::bad_alloc &) {
		releaseResources();

		return RetCode::E_NO_MEMORY;
	}
	ionicStrength = internal->rawIonicStrength();
	if (iterationsNeeded != nullptr)
		*iterationsNeeded = internal->iterations();

	releaseResources();

	return RetCode::OK;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_SOLVERIMPL_HPP
