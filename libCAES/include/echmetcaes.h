#ifndef ECHMET_CAES_CAES_H
#define ECHMET_CAES_CAES_H

#include <cstddef>

#ifdef ECHMET_HAVE_CPP11
#include <cstdint>
#else
#include <stdint.h>
#endif // ECHMET_HAVE_CPP11

#define ECHMET_IMPORT_INTERNAL
#include <echmetsyscomp.h>
#undef ECHMET_IMPORT_INTERNAL

#include <echmetmodule.h>

namespace ECHMET {

/*!
 * Complexation and Acidobazic Equilibria Solver package.
 */
namespace CAES {

/*!
 * Number of iterations done by the numerical solver.
 */
class SolverIterations {
public:
	uint32_t outer;	/*!< Number of outer loop iterations. The outer loop runs multiple times only when correction for ionic strength is requested. */
	uint32_t total; /*!< Total number of times the inner loop was run. */

};
IS_POD(SolverIterations)

/*!
 * SolverContext interface.
 */
class SolverContext {
public:
	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~SolverContext() ECHMET_NOEXCEPT = 0;
};

/*!
 * Equilibria solver interface.
 */
class Solver {
public:
	/*!
	 * Options that modify behavior of the solver.
	 */
	ECHMET_WK_ENUM(Options) {
		NONE = 0,
		DISABLE_THREAD_SAFETY = (1 << 0)	/*!< Use thread-unsafe variant of the solver.
							     This may improve performance in cases where the solver is not
							     used from multiple threads.
							     This option cannot be changed during solver's lifetime. */
		ENUM_FORCE_INT32_SIZE(CAESOptions)
	};

	/*!
	 * Returns a pointer to the context assigned to the solver.
	 *
	 * @return Pointer to \p SolverContext object.
	 */
	virtual SolverContext * ECHMET_CC context() ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns default solver options.
	 *
	 * @return Default options.
	 */
	static Options ECHMET_CC defaultOptions() ECHMET_NOEXCEPT
	{
		return static_cast<Options>(0);
	}

	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Calculates the initial estimation of concentration of all species in the system.
	 * This is a fast estimation function that requires a "reasonably good" initial estimate
	 * of the \p cH concentration. The function might return \p RetCode::OK even if the calculated
	 * distribution is incorrect. Use this carefully.
	 *
	 * @param[in] cHInitial Initial estimate of H<sub>3</sub>O<sup>+</sup> ions.
	 * @param[in] analyticalConcentrations Vector of analytical concentrations of all compounds in the system.
	 * @param[in,out] calcProps \p CalculatedProperties object associated with the system that is being solved.
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of concentrations vector.
	 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution.
	 * @retval RetCode::E_FAST_ESTIMATE_FAILURE Fast estimation failed to find a solution.
	 */
	virtual RetCode ECHMET_CC estimateDistributionFast(const ECHMETReal &cHInitial, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) ECHMET_NOEXCEPT = 0;

	/*!
	 * Calculates the initial estimation of concentration of all species in the system.
	 * This is a safe estimation function that is guaranteed to converge for any range of pH
	 * between -3 to \p inf
	 *
	 * @param[in] analyticalConcentrations Vector of analytical concentrations of all compounds in the system.
	 * @param[in,out] calcProps \p CalculatedProperties object associated with the system that is being solved.
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of concentrations vector.
	 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution.
	 */
	virtual RetCode ECHMET_CC estimateDistributionSafe(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns the options currently set for the solver.
	 *
	 * @return \p Options enum.
	 */
	virtual Options ECHMET_CC options() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Sets a new context for the solver.
	 *
	 * @param[in] ctx New context.
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_INVALID_ARGUMENT \p SolverContext not castable to \p SolverContextImpl.
	 * @retval RetCode::E_NO_MEMORY Not enough memory to allocate new \p SolverInternal object.
	 */
	virtual RetCode ECHMET_CC setContext(SolverContext *ctx) ECHMET_NOEXCEPT = 0;

	/*!
	 * Sets new options for the solver.
	 *
	 * @param[in] options The new options to set.
	 *
	 * @retval RetCode::OK SUCCESS
	 * @retval RetCode::E_INVALID_ARGUMENT Nonsensical options passed as the argument
	 */
	virtual RetCode ECHMET_CC setOptions(const Options options) ECHMET_NOEXCEPT = 0;

	/*!
	 * Calculates the equilibrium ionic distribution of the system.
	 * The ionic composition shall be estimate by calling \p estimateDistribution() prior to calling \p solve().
	 *
	 * @param[in] analyticalConcentrations Analytical concentrations of constituents in the system.
	 * @param[in,out] calcProps Struct where to store the results.
	 * @param[in] iterations Maximum number of iterations to try.
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
	 * @retval RetCode::E_INVALID_ARGUMENT Initial value of ionic strength is invalid.
	 */
	virtual RetCode ECHMET_CC solve(const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const size_t iterations, SolverIterations *iterationsNeeded = ECHMET_NULLPTR) ECHMET_NOEXCEPT = 0;

protected:
	virtual ~Solver() ECHMET_NOEXCEPT = 0;
};

extern "C" {

/*!
 * Create instance of CAES Solver with double precision.
 *
 * @param[in] ctx Contexx for the solver. The context must be of the double-precision type
 * @param[in] options Solver options
 * @param[in] corrections Corrections for non-ideal solutions to use
 *
 * @return Pointer to the Solver instance or \p NULL if the Solver could not have been created
 */
ECHMET_API Solver * ECHMET_CC createSolver(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) ECHMET_NOEXCEPT;

/*!
 * Create instance of CAES Solver with high precision.
 *
 * @param[in] ctx Contexx for the solver. The context must be of the high-precision type
 * @param[in] options Solver options
 * @param[in] corrections Corrections for non-ideal solutions to use
 *
 * @return Pointer to the Solver instance or \p NULL if the Solver could not have been created
 */
ECHMET_API Solver * ECHMET_CC createSolverHighPrecision(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) ECHMET_NOEXCEPT;

/*!
 * Creates a solver context for a given system and intializes
 * total and ionic concentrations vectors and their mappings.
 * This creates context for double precision solver.
 *
 * @param[in,out] ctx Reference to SolverContext to be created by this function.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Not enough memory to create the sovler context.
 * @retval RetCode::E_DATA_TOO_LARGE Amount of data to be processed is too large.
 * @retval RetCode::E_BAD_INPUT Nonsensical input data.
 * @retval RetCode::E_MISSING_PB Complexation constant was not set.
 */
ECHMET_API RetCode ECHMET_CC createSolverContext(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) ECHMET_NOEXCEPT;

/*!
 * Creates a solver context for a given system and intializes
 * total and ionic concentrations vectors and their mappings.
 * This creates context for high precision solver.
 *
 * @param[in,out] ctx Reference to SolverContext to be created by this function.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Not enough memory to create the sovler context.
 * @retval RetCode::E_DATA_TOO_LARGE Amount of data to be processed is too large.
 * @retval RetCode::E_BAD_INPUT Nonsensical input data.
 * @retval RetCode::E_MISSING_PB Complexation constant was not set.
 */
ECHMET_API RetCode ECHMET_CC createSolverContextHighPrecision(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) ECHMET_NOEXCEPT;

} // extern "C"

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CAES_H
