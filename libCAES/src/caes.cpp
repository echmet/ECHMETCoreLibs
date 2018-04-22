#include "caes_p.h"

namespace ECHMET {
namespace CAES {

Solver * ECHMET_CC createSolver(SolverContext *ctx, const NonidealityCorrections corrections) noexcept
{
	return createSolverInternal<ECHMETReal>(ctx, corrections);
}

/*!
 * Creates a solver context for a given system and intializes
 * total and ionic concentrations vectors and their mappings.
 *
 * @param[in,out] ctx Reference to SolverContext to be created by this function.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Not enough memory to create the sovler context.
 * @retval RetCode::E_DATA_TOO_LARGE Amount of data to be processed is too large.
 * @retval RetCode::E_BAD_INPUT Nonsensical input data.
 * @retval RetCode::E_MISSING_PB Complexation constant was not set.
 */
RetCode ECHMET_CC createSolverContext(SolverContext *&ctx, const SolverContext::Options options, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	return createSolverContextInternal<ECHMETReal>(ctx, options, chemSystem);
}

/*!
 * Calculates the initial estimation of concentration of all species in the system
 *
 * @param[in] ctx Solver context
 * @param[in,out] ionicConcentrations Vector of concentrations of all ionic forms
 *                The vector shall have the expected size
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_INVALID_ARGUMENT Unexpected size of the concentration vectors or the \ctx
 *         pointer is not castable to \SolverContextImpl
 * @retval RetCode::E_NO_MEMORY Not enough memory to estimate distribution
 */
RetCode estimateDistribution(SolverContext *ctx, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	SolverVector<ECHMETReal> estC{};

	RetCode tRet = estimateDistributionInternal<ECHMETReal>(ctx, analyticalConcentrations, estC);
	if (tRet != RetCode::OK)
		return tRet;

	for (int idx = 0; idx < estC.size(); idx++)
		(*calcProps.ionicConcentrations)[idx] = estC(idx);

	return RetCode::OK;
}

SolverContext::~SolverContext() noexcept {}
Solver::~Solver() noexcept {}

} // namespace CAES

} // namespace ECHMET
