#ifndef ECHMET_CAES_SOLVERIMPL_HPP
#define ECHMET_CAES_SOLVERIMPL_HPP

#include "caes_p.h"
#include "solverinternal.h"
#include "estimator_helpers.hpp"

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

InstructionSet detectInstructionSet();

template <typename CAESReal, bool V = std::is_same<CAESReal, mpfr::mpreal>::value>
static
FreeMPFRCacheSwitch<V> freeMPFRCache()
{
	return FreeMPFRCacheSwitch<V>{};
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept
{
	if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY)) {
		delete internal;
		mpfr_free_cache();
	}
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept
{
	if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY))
		delete internal;
}

/*!
 * Solver c-tor.
 *
 * @param[in] ctx \p SolverContext used to initialize solver's internal data structures. The pointer
 *            shall remain valid throughout the entire lifetime of the \p Solver .
 * @param[in] options Solver options.
 */
template <typename CAESReal>
SolverImpl<CAESReal>::SolverImpl(SolverContextImpl<CAESReal> *ctx, const Options options, const NonidealityCorrections corrections) :
	m_ctx(ctx),
	m_options(options),
	m_internalUnsafe(nullptr),
	m_anCVecUnsafe(nullptr),
	m_estimatedConcentrationsUnsafe(nullptr),
	m_instructionSet(detectInstructionSet())
{
	m_correctDebyeHuckel = nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

	initializeTotalEquilibria(ctx);

	if (m_options & Solver::Options::DISABLE_THREAD_SAFETY) {
		try {
			initializeEstimators();

			m_internalUnsafe = makeSolverInternal(ctx);
			m_anCVecUnsafe = new SolverVector<CAESReal>(ctx->analyticalConcentrationCount);
			m_estimatedConcentrationsUnsafe = AlignedAllocator<CAESReal, 32>::alloc(ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete m_internalUnsafe;
			delete m_anCVecUnsafe;
			releaseTotalEquilibria();

			throw;
		}
	}
}

/*!
 * Solver d-tor.
 */
template <typename CAESReal>
SolverImpl<CAESReal>::~SolverImpl() ECHMET_NOEXCEPT
{
	delete m_internalUnsafe;
	delete m_anCVecUnsafe;

	releaseTotalEquilibria();
	AlignedAllocator<CAESReal, 32>::free(m_estimatedConcentrationsUnsafe);
}

/*!
 * Returns a pointer to \p SolverContext assigned to the solve.
 *
 * @return Pointer to \p SolverContext object.
 */
template <typename CAESReal>
SolverContext * ECHMET_CC SolverImpl<CAESReal>::context() ECHMET_NOEXCEPT
{
	return m_ctx;
}

template <typename CAESReal>
void ECHMET_CC SolverImpl<CAESReal>::destroy() const ECHMET_NOEXCEPT
{
	delete this;
}

template <typename CAESReal>
void SolverImpl<CAESReal>::defaultActivityCoefficients(std::vector<CAESReal> &activityCoefficients) const
{
	std::fill(activityCoefficients.begin(), activityCoefficients.end(), 1.0);
}

template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::estimateDistributionFast(const ECHMETReal &cHInitial, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	CAESReal ionicStrength;
	SolverVector<CAESReal> estC{};

	RetCode tRet = estimateDistributionInternal(cHInitial, analyticalConcentrations, estC, ionicStrength, true);
	if (tRet != RetCode::OK)
		return tRet;

	for (int idx = 0; idx < estC.size(); idx++)
		(*calcProps.ionicConcentrations)[idx] = CAESRealToECHMETReal(estC(idx));
	calcProps.ionicStrength = CAESRealToECHMETReal(ionicStrength);

	return RetCode::OK;
}

template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::estimateDistributionSafe(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	const ECHMETReal dummy = 0.0;	/* Unused */
	CAESReal ionicStrength;
	SolverVector<CAESReal> estC{};

	RetCode tRet = estimateDistributionInternal(dummy, analyticalConcentrations, estC, ionicStrength, false);
	if (tRet != RetCode::OK)
		return tRet;

	for (int idx = 0; idx < estC.size(); idx++)
		(*calcProps.ionicConcentrations)[idx] = CAESRealToECHMETReal(estC(idx));
	calcProps.ionicStrength = CAESRealToECHMETReal(ionicStrength);

	return RetCode::OK;
}

/*!
 * Calculates the initial estimation of concentration of all species in the system.
 *
 * @param[in] ctx Solver context.
 * @param[in,out] ionicConcentrations Vector of concentrations of all ionic forms.
 *                The vector shall have the expected size.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of the concentration vectors.
 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution.
 * @retval RetCode::E_FAST_ESTIMATE_FAILURE Fast estimation failed to find a solution.
 */
template <typename CAESReal>
RetCode SolverImpl<CAESReal>::estimateDistributionInternal(const CAESReal &cHInitial, const RealVec *analyticalConcentrations, SolverVector<CAESReal> &estimatedConcentrations,
							   CAESReal &ionicStrength, const bool useFastEstimate) noexcept
{
	estimatedConcentrations.resize(m_ctx->concentrationCount);

	try {
		const auto results = [&]() {
			if (useFastEstimate) {
				if (m_options & Solver::Options::DISABLE_THREAD_SAFETY) {
					return estimatepHFast<false>(cHInitial, analyticalConcentrations,
								     m_estimatedIonicConcentrations, m_dEstimatedIonicConcentrationsdH,
								     m_activityCoefficients);
				} else {
					std::vector<CAESReal> activityCoefficients;
					SolverVector<CAESReal> icConcs(m_TECount);
					SolverVector<CAESReal> dIcConcsdH(m_TECount);

					activityCoefficients.resize(m_ctx->chargesSquared.size());
					return estimatepHFast<true>(cHInitial, analyticalConcentrations, icConcs, dIcConcsdH, m_activityCoefficients);
				}
			} else {
				if (m_options & Solver::Options::DISABLE_THREAD_SAFETY)
					return estimatepHSafe<false>(analyticalConcentrations, m_estimatedIonicConcentrations, m_activityCoefficients);
				else {
					std::vector<CAESReal> activityCoefficients;
					SolverVector<CAESReal> icConcs(m_TECount);

					activityCoefficients.resize(m_ctx->chargesSquared.size());
					return estimatepHSafe<true>(analyticalConcentrations, icConcs, activityCoefficients);
				}
			}
		}();
		ionicStrength = results.second;

		ECHMET_DEBUG_CODE(for (int idx = 0; idx < results.first.size(); idx++) {
				  const CAESReal &v = results.first(idx);
				  fprintf(stderr, "estC: %.4g, pX: %.4g\n", CAESRealToDouble(v), CAESRealToDouble(pX(v)));
				  });

		ECHMET_DEBUG_CODE(fprintf(stderr, "Estimated pH = %g\n", CAESRealToDouble(pX(results.first(0)) + 3.0)));
		ECHMET_DEBUG_CODE(fprintf(stderr, "Estimated ionic strength = %g\n", CAESRealToDouble(ionicStrength)));
		/* H+ and OH- are expected to be the first and second item in the vector */
		estimatedConcentrations(0) = results.first(0);
		estimatedConcentrations(1) = results.first(1);

		estimateComplexesDistribution<CAESReal>(m_ctx->complexNuclei, m_ctx->allLigands, results.first, m_ctx->allForms->size() + 2, estimatedConcentrations);
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
template <typename CAESReal> template <bool ThreadSafe>
std::pair<SolverVector<CAESReal>, CAESReal> SolverImpl<CAESReal>::estimatepHFast(const CAESReal &cHInitial, const RealVec *analyticalConcentrations,
										 SolverVector<CAESReal> &icConcs, SolverVector<CAESReal> &dIcConcsdH,
										 std::vector<CAESReal> &activityCoefficients)
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

		while (true) {
			calculateDistributionWithDerivative<CAESReal, ThreadSafe>(cH, icConcs, dIcConcsdH, m_totalEquilibria, analyticalConcentrations, activityCoefficients);

			CAESReal z = calcTotalCharge<CAESReal, ThreadSafe>(icConcs, m_totalEquilibria);
			CAESReal dZ = calcTotalCharge<CAESReal, ThreadSafe>(dIcConcsdH, m_totalEquilibria);

			cOH = KW_298 / (cH * activityOneSquared);

			z += cH - cOH;
			dZ += 1.0 + cOH * cOH;

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
		icConcs(1) = cOH;

		ionicStrength = calculateIonicStrength<CAESReal, ThreadSafe>(icConcs, m_totalEquilibria, m_ctx->chargesSquared);
		if (m_correctDebyeHuckel) {
			calculateActivityCoefficients(ionicStrength, activityCoefficients, m_ctx->chargesSquared);
			ionicStrength = !hasIonicStrengthConverged(maxChargeActivityCoeff, activityCoefficients.back());
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
template <typename CAESReal> template <bool ThreadSafe>
std::pair<SolverVector<CAESReal>, CAESReal> SolverImpl<CAESReal>::estimatepHSafe(const RealVec *analyticalConcentrations, SolverVector<CAESReal> &icConcs, std::vector<CAESReal> &activityCoefficients)
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
			calculateDistribution<CAESReal, ThreadSafe>(cH, icConcs, m_totalEquilibria, analyticalConcentrations, activityCoefficients);

			CAESReal z = calcTotalCharge<CAESReal, ThreadSafe>(icConcs, m_totalEquilibria);

			cOH = KW_298 / (cH * activityOneSquared);
			z += cH - cOH;
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

		icConcs(0) = cH;
		icConcs(1) = cOH;

		ionicStrength = calculateIonicStrength<CAESReal, ThreadSafe>(icConcs, m_totalEquilibria, m_ctx->chargesSquared);
		if (m_correctDebyeHuckel) {
			calculateActivityCoefficients(ionicStrength, activityCoefficients, m_ctx->chargesSquared);
			ionicStrength = !hasIonicStrengthConverged(maxChargeActivityCoeff, activityCoefficients.back());
		} else
			ionicStrengthUnstable = false;
	} while (m_correctDebyeHuckel && ionicStrengthUnstable && isLoopCtr++ < 100);

	return { icConcs, ionicStrength };
}

template <typename CAESReal>
void SolverImpl<CAESReal>::initializeEstimators()
{
	/* The +2 is for H+ and OH- */
	m_estimatedIonicConcentrations.resize(m_TECount);
	m_dEstimatedIonicConcentrationsdH.resize(m_TECount);
	m_activityCoefficients.resize(m_ctx->chargesSquared.size());
	defaultActivityCoefficients(m_activityCoefficients);
}

/*!
 * Initializes TotalEquilibria objects used to estimate ionic distribution
 *
 * @param[in] ctx SovlerContextImpl to use
 */
template <typename CAESReal>
void SolverImpl<CAESReal>::initializeTotalEquilibria(const SolverContextImpl<CAESReal> *ctx)
{
	auto countTEForms = [](const auto &TEVec, const auto &caster) {
		size_t TECount = 0;
		for (const auto *teb : TEVec) {
			const auto *te = caster(teb);
			TECount += te->numHigh - te->numLow + 1;
		}
		return TECount;
	};

	auto makeTotalEquilibrium = [this](const auto *ct) -> TotalEquilibriumBase * {
		if (m_options & Solver::Options::DISABLE_THREAD_SAFETY)
			return new TotalEquilibrium<CAESReal, false>(ct->chargeLow, ct->chargeHigh, ct->pKas, ct->analyticalConcentrationIndex);
		return new TotalEquilibrium<CAESReal, true>(ct->chargeLow, ct->chargeHigh, ct->pKas, ct->analyticalConcentrationIndex);
	};

	releaseTotalEquilibria();
	m_totalEquilibria.reserve(ctx->complexNuclei->size() + ctx->allLigands->size());

	for (const ComplexNucleus<CAESReal> *cn : *ctx->complexNuclei)
		m_totalEquilibria.emplace_back(makeTotalEquilibrium(cn));
	for (const Ligand<CAESReal> *l : *ctx->allLigands)
		m_totalEquilibria.emplace_back(makeTotalEquilibrium(l));

	if (m_options & Solver::Options::DISABLE_THREAD_SAFETY)
		m_TECount = countTEForms(m_totalEquilibria, [](const TotalEquilibriumBase *teb) { return static_cast<const TotalEquilibrium<CAESReal, false> *>(teb); });
	else
		m_TECount = countTEForms(m_totalEquilibria, [](const TotalEquilibriumBase *teb) { return static_cast<const TotalEquilibrium<CAESReal, true> *>(teb); });
	m_TECount += 2;
}

/*!
 * Makes appropriate \p SolverInternal object
 * based on available SIMD instruction set
 *
 * @retval Pointer to \p SolverInternal object
 */
template <typename CAESReal>
SolverInternalBase<CAESReal> * SolverImpl<CAESReal>::makeSolverInternal(const SolverContextImpl<CAESReal> *ctx) const
{
	switch (m_instructionSet) {
	case InstructionSet::SSE2:
		return new SolverInternal<CAESReal, InstructionSet::SSE2>(ctx);
	case InstructionSet::AVX:
		return new SolverInternal<CAESReal, InstructionSet::AVX>(ctx);
	case InstructionSet::FMA3:
		return new SolverInternal<CAESReal, InstructionSet::FMA3>(ctx);
	default:
		return new SolverInternal<CAESReal, InstructionSet::GENERIC>(ctx);
	}
}

/*!
 * Returns the current options of the solver.
 *
 * @return Solver options.
 */
template <typename CAESReal>
typename SolverImpl<CAESReal>::Options ECHMET_CC SolverImpl<CAESReal>::options() const noexcept
{
	return m_options;
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseTotalEquilibria()
{
	for (auto teb : m_totalEquilibria)
		delete teb;

	m_totalEquilibria.clear();
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
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::setContext(SolverContext *ctx) ECHMET_NOEXCEPT
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
template <typename CAESReal>
RetCode SolverImpl<CAESReal>::setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept
{
	releaseTotalEquilibria();

	try {
		initializeTotalEquilibria(ctx);
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	if (m_options & Solver::Options::DISABLE_THREAD_SAFETY) {
		delete m_internalUnsafe;
		delete m_anCVecUnsafe;
		AlignedAllocator<CAESReal, 32>::free(m_estimatedConcentrationsUnsafe);

		m_internalUnsafe = nullptr;
		m_anCVecUnsafe = nullptr;
		m_estimatedConcentrationsUnsafe = nullptr;

		try {
			initializeEstimators();
			m_internalUnsafe = makeSolverInternal(ctx);
			m_anCVecUnsafe = new SolverVector<CAESReal>(ctx->analyticalConcentrationCount);
			m_estimatedConcentrationsUnsafe = AlignedAllocator<CAESReal, 32>::alloc(ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete m_anCVecUnsafe;
			delete m_internalUnsafe;
			releaseTotalEquilibria();

			return RetCode::E_NO_MEMORY;
		}
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
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::setOptions(const Options options) noexcept
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
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded) ECHMET_NOEXCEPT
{
	SolverInternalBase<CAESReal> *internal = nullptr;
	SolverVector<CAESReal> *anCVec = nullptr;
	CAESReal *estimatedConcentrations = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY)) {
			delete anCVec;
			AlignedAllocator<CAESReal, 32>::free(estimatedConcentrations);
		}
	};

	if (m_options & Solver::Options::DISABLE_THREAD_SAFETY) {
		internal = m_internalUnsafe;
		anCVec = m_anCVecUnsafe;
		estimatedConcentrations = m_estimatedConcentrationsUnsafe;
	} else {
		try {
			internal = makeSolverInternal(m_ctx);
			anCVec = new SolverVector<CAESReal>(m_ctx->analyticalConcentrationCount);
			estimatedConcentrations = AlignedAllocator<CAESReal, 32>::alloc(m_ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete anCVec;
			delete internal;
			return RetCode::E_NO_MEMORY;
		}
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
template <typename CAESReal>
RetCode SolverImpl<CAESReal>::solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const SolverVector<CAESReal> &estimatedConcentrations, const size_t iterations, SolverIterations *iterationsNeeded) noexcept
{
	SolverInternalBase<CAESReal> *internal = nullptr;
	CAESReal *estimatedConcentrationsInternal = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (!(m_options & Solver::Options::DISABLE_THREAD_SAFETY))
			AlignedAllocator<CAESReal, 32>::free(estimatedConcentrationsInternal);
	};

	if (m_options & Solver::Options::DISABLE_THREAD_SAFETY) {
		internal = m_internalUnsafe;
		estimatedConcentrationsInternal = m_estimatedConcentrationsUnsafe;
	} else {
		try {
			internal = makeSolverInternal(m_ctx);
			estimatedConcentrationsInternal = AlignedAllocator<CAESReal, 32>::alloc(m_ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete internal;
			return RetCode::E_NO_MEMORY;
		}
	}

	for (int idx = 0; idx < estimatedConcentrations.rows(); idx++)
		estimatedConcentrationsInternal[idx] = estimatedConcentrations(idx);

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
