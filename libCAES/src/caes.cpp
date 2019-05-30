#include "caes_p.h"

#define ECHMET_IMPORT_INTERNAL
#include <echmetelems.h>

#ifdef ECHMET_DEBUG_OUTPUT
#include <cstdio>
#endif // ECHMET_DEBUG_OUTPUT

namespace ECHMET {
namespace CAES {

Solver * ECHMET_CC createSolver(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) noexcept
{
	return createSolverInternal<ECHMETReal>(ctx, options, corrections);
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
RetCode ECHMET_CC createSolverContext(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	return createSolverContextInternal<ECHMETReal>(ctx, chemSystem);
}

InstructionSet detectInstructionSet() noexcept
{
	const CPUSIMD simd = cpuSupportedSIMD();

	if (simd.FMA3) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Using FMA3-optimized solver\n"));
		return InstructionSet::FMA3;
	}
	if (simd.AVX) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Using AVX-optimized solver\n"));
		return InstructionSet::AVX;
	}
	if (simd.SSE2) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Using SSE2-optimized solver\n"));
		return InstructionSet::SSE2;
	}

	ECHMET_DEBUG_CODE(fprintf(stderr, "Using generic solver\n"));

	return InstructionSet::GENERIC;
}

SolverContext::~SolverContext() noexcept {}
Solver::~Solver() noexcept {}

} // namespace CAES

} // namespace ECHMET
