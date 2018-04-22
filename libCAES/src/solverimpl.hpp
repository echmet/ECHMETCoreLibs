#ifndef ECHMET_CAES_SOLVERIMPL_HPP
#define ECHMET_CAES_SOLVERIMPL_HPP

#include "caes_p.h"
#include "solverinternal.h"

#include <x86intrin.h>

namespace ECHMET {
namespace CAES {

InstructionSet detectInstructionSet();

template <typename CAESReal, bool V = std::is_same<CAESReal, mpfr::mpreal>::value>
FreeMPFRCacheSwitch<V> freeMPFRCache()
{
	return FreeMPFRCacheSwitch<V>{};
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept
{
	if (!(m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY)) {
		delete internal;
		mpfr_free_cache();
	}
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternalBase<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept
{
	if (!(m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY))
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
SolverImpl<CAESReal>::SolverImpl(SolverContextImpl<CAESReal> *ctx, const NonidealityCorrections corrections) :
	m_ctx(ctx),
	m_internalUnsafe(nullptr),
	m_anCVecUnsafe(nullptr),
	m_estimatedConcentrationsUnsafe(nullptr),
	m_instructionSet(detectInstructionSet())
{
	m_correctDebyeHuckel = nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

	if (m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) {
		try {
			m_anCVecUnsafe = new SolverVector<CAESReal>(ctx->analyticalConcentrationCount);
			m_estimatedConcentrationsUnsafe = alignedAlloc<CAESReal>(ctx->concentrationCount);
			m_internalUnsafe = makeSolverInternal(ctx);
		} catch (const std::bad_alloc &up) {
			delete m_anCVecUnsafe;
			alignedFree(m_estimatedConcentrationsUnsafe);

			throw up;
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
	alignedFree(m_estimatedConcentrationsUnsafe);
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
 * Sets new \p SolverContext.
 *
 * @param[in] ctx Pointer to the new \p SolverContext.
 * @param[in] totalConcentrations Pointer to vector of new total concentrations of all system constituents.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT SolverContext not castable to SolverContextImpl.
 * @retval RetCode::E_INVALID_ARGUMENT Change of thread safetiness is not allowed.
 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate new \p SolverInternal object.
 */
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::setContext(SolverContext *ctx) ECHMET_NOEXCEPT
{
	SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	if ((m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) !=
	    (ctxImpl->options() & SolverContext::Options::DISABLE_THREAD_SAFETY))
		return RetCode::E_INVALID_ARGUMENT;

	return setContextInternal(ctxImpl);
}

/*!
 * Sets new \p SolverContext. Internal imlementation.
 *
 * @param[in] ctx Pointer to the new \p SolverContext.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate new \p SolverInternal object.
 */
template <typename CAESReal>
RetCode SolverImpl<CAESReal>::setContextInternal(SolverContextImpl<CAESReal> *ctx) noexcept
{
	if (ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) {
		delete m_internalUnsafe;
		delete m_anCVecUnsafe;
		alignedFree(m_estimatedConcentrationsUnsafe);

		m_internalUnsafe = nullptr;
		m_anCVecUnsafe = nullptr;
		m_estimatedConcentrationsUnsafe = nullptr;

		try {
			m_internalUnsafe = makeSolverInternal(ctx);
			m_anCVecUnsafe = new SolverVector<CAESReal>(ctx->analyticalConcentrationCount);
			m_estimatedConcentrationsUnsafe = alignedAlloc<CAESReal>(ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete m_anCVecUnsafe;
			delete m_internalUnsafe;

			return RetCode::E_NO_MEMORY;
		}
	}
	m_ctx = ctx;

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
 * @retval RetCode::E_INVALID_ARGUMENT Vector of analytical concentrations has invalid size.
 */
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded) ECHMET_NOEXCEPT
{
	SolverInternalBase<CAESReal> *internal = nullptr;
	SolverVector<CAESReal> *anCVec = nullptr;
	CAESReal *estimatedConcentrations = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (!(m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY)) {
			delete anCVec;
			alignedFree(estimatedConcentrations);
		}
	};

	if (m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) {
		internal = m_internalUnsafe;
		anCVec = m_anCVecUnsafe;
		estimatedConcentrations = m_estimatedConcentrationsUnsafe;
	} else {
		try {
			internal = makeSolverInternal(m_ctx);
			anCVec = new SolverVector<CAESReal>(m_ctx->analyticalConcentrationCount);
			estimatedConcentrations = alignedAlloc<CAESReal>(m_ctx->concentrationCount);
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

		const RetCode tRet = internal->solve(anCVec, estimatedConcentrations, m_correctDebyeHuckel, iterations);
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
 * @param[out] ionicStrength Ionic strength of the solved system.
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
 */
template <typename CAESReal>
RetCode SolverImpl<CAESReal>::solveRaw(SolverVector<CAESReal> &concentrations, CAESReal &ionicStrength, const SolverVector<CAESReal> *anCVec, const SolverVector<CAESReal> &estimatedConcentrations, const size_t iterations, SolverIterations *iterationsNeeded) noexcept
{
	SolverInternalBase<CAESReal> *internal = nullptr;
	CAESReal *estimatedConcentrationsInternal = nullptr;

	const auto releaseResources = [&]() {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());
		if (!(m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY))
			alignedFree(estimatedConcentrationsInternal);
	};

	if (m_ctx->options() & SolverContext::Options::DISABLE_THREAD_SAFETY) {
		internal = m_internalUnsafe;
		estimatedConcentrationsInternal = m_estimatedConcentrationsUnsafe;
	} else {
		try {
			internal = makeSolverInternal(m_ctx);
			estimatedConcentrationsInternal = alignedAlloc<CAESReal>(m_ctx->concentrationCount);
		} catch (const std::bad_alloc &) {
			delete internal;
			return RetCode::E_NO_MEMORY;
		}
	}

	for (int idx = 0; idx < estimatedConcentrations.rows(); idx++)
		estimatedConcentrationsInternal[idx] = estimatedConcentrations(idx);

	const RetCode tRet = internal->solve(anCVec, estimatedConcentrationsInternal, m_correctDebyeHuckel, iterations);
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
