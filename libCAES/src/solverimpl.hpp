#ifndef ECHMET_CAES_SOLVERIMPL_HPP
#define ECHMET_CAES_SOLVERIMPL_HPP

#include "caes_p.h"
#include "solverinternal.h"

namespace ECHMET {
namespace CAES {

template <typename CAESReal, bool V = std::is_same<CAESReal, mpfr::mpreal>::value>
FreeMPFRCacheSwitch<V> freeMPFRCache()
{
	return FreeMPFRCacheSwitch<V>{};
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternal<CAESReal> *internal, FreeMPFRCacheSwitch<true>) noexcept
{
	if (!(m_options & Options::DISABLE_THREAD_SAFETY)) {
		delete internal;
		mpfr_free_cache();
	}
}

template <typename CAESReal>
void SolverImpl<CAESReal>::releaseSolverInternal(SolverInternal<CAESReal> *internal, FreeMPFRCacheSwitch<false>) noexcept
{
	if (!(m_options & Options::DISABLE_THREAD_SAFETY))
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
SolverImpl<CAESReal>::SolverImpl(const SolverContextImpl<CAESReal> *ctx, const NonidealityCorrections corrections, const Options options) noexcept :
	m_options(options),
	m_ctx(ctx),
	m_internalUnsafe(nullptr)
{
	m_correctDebyeHuckel = nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

	if (m_options & Options::DISABLE_THREAD_SAFETY)
		m_internalUnsafe = new (std::nothrow) SolverInternal<CAESReal>(ctx);
}

/*!
 * Solver d-tor.
 */
template <typename CAESReal>
SolverImpl<CAESReal>::~SolverImpl() ECHMET_NOEXCEPT
{
	delete m_internalUnsafe;
}

/*!
 * Returns a pointer to \p SolverContext assigned to the solve.
 *
 * @return Pointer to \p SolverContext object.
 */
template <typename CAESReal>
const SolverContext * ECHMET_CC SolverImpl<CAESReal>::context() const ECHMET_NOEXCEPT
{
	return m_ctx;
}

template <typename CAESReal>
void ECHMET_CC SolverImpl<CAESReal>::destroy() const ECHMET_NOEXCEPT
{
	delete this;
}

/*!
 * Returns the current options of the solver.
 *
 * @return Solver options.
 */
template <typename CAESReal>
typename SolverImpl<CAESReal>::Options ECHMET_CC SolverImpl<CAESReal>::options() const ECHMET_NOEXCEPT
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
 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate new \p SolverInternal object.
 */
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::setContext(const SolverContext *ctx) ECHMET_NOEXCEPT
{
	const SolverContextImpl<CAESReal> *ctxImpl = dynamic_cast<const SolverContextImpl<CAESReal> *>(ctx);
	if (ctxImpl == nullptr)
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
RetCode SolverImpl<CAESReal>::setContextInternal(const SolverContextImpl<CAESReal> *ctx) noexcept
{
	if (m_options & Options::DISABLE_THREAD_SAFETY) {
		delete m_internalUnsafe;

		m_internalUnsafe = new (std::nothrow) SolverInternal<CAESReal>(ctx);

		if (m_internalUnsafe == nullptr)
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
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::setOptions(const Options options) ECHMET_NOEXCEPT
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
 * @retval RetCode::E_SOLVER_NOT_INITIALIZED Solver was not initialized prior to calling this function.
 * @retval RetCode::E_NO_MEMORY Not enough memory to perform the calculation.
 * @retval RetCode::E_NRS_FAILURE Newton-Rapshon solver encountered an error during calculation.
 * @retval RetCode::E_NRS_NO_CONVERGENCE Newton-Raphson solver failed to converge within the given number of iterations
 * @retval RetCode::E_NRS_STUCK Greatest change of X-value calculated by the Newton-Raphson solver is below the precision threshold.
 * @retval RetCode::E_NRS_NO_SOLUTION System appears to have no solution.
 * @retval RetCode::E_IS_NO_CONVERGENCE Solver failed to find a solution within the given number of iterations.
 */
template <typename CAESReal>
RetCode ECHMET_CC SolverImpl<CAESReal>::solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded) ECHMET_NOEXCEPT
{
	SolverInternal<CAESReal> *internal = nullptr;
	if (m_options & Options::DISABLE_THREAD_SAFETY)
		internal = m_internalUnsafe;
	else
		internal = new (std::nothrow) SolverInternal<CAESReal>(m_ctx);

	if (internal == nullptr)
		return RetCode::E_SOLVER_NOT_INITIALIZED;

	try {
		SolverVector<CAESReal> anCVec{analyticalConcentrations->size()};
		SolverVector<CAESReal> estimatedConcentrations{calcProps.ionicConcentrations->size()};

		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			anCVec(idx) = analyticalConcentrations->at(idx);

		for (size_t idx = 0; idx < calcProps.ionicConcentrations->size(); idx++)
			estimatedConcentrations(idx) = calcProps.ionicConcentrations->at(idx);

		const RetCode tRet = internal->solve(&anCVec, estimatedConcentrations, m_correctDebyeHuckel, iterations);
		if (tRet != RetCode::OK) {
			releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

			return tRet;
		}

		internal->resultsToOutput(calcProps);

		if (iterationsNeeded != nullptr)
			*iterationsNeeded = internal->iterations();

		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

		return RetCode::OK;
	} catch (std::bad_alloc &) {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

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
	SolverInternal<CAESReal> *internal = nullptr;
	if (m_options & Options::DISABLE_THREAD_SAFETY)
		internal = m_internalUnsafe;
	else
		internal = new (std::nothrow) SolverInternal<CAESReal>(m_ctx);

	if (internal == nullptr)
		return RetCode::E_SOLVER_NOT_INITIALIZED;

	const RetCode tRet = internal->solve(anCVec, estimatedConcentrations, m_correctDebyeHuckel, iterations);
	if (tRet != RetCode::OK) {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

		return tRet;
	}

	try {
		concentrations = internal->rawConcentrations();
	} catch (std::bad_alloc &) {
		releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

		return RetCode::E_NO_MEMORY;
	}
	ionicStrength = internal->rawIonicStrength();
	if (iterationsNeeded != nullptr)
		*iterationsNeeded = internal->iterations();

	releaseSolverInternal(internal, freeMPFRCache<CAESReal>());

	return RetCode::OK;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_SOLVERIMPL_HPP
