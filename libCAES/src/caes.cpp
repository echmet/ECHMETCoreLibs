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

Solver * ECHMET_CC createSolverHighPrecision(SolverContext *ctx, const Solver::Options options, const NonidealityCorrections corrections) noexcept
{
	NumericPrecisionSetter<mpfr::mpreal> nps{};

	return createSolverInternal<mpfr::mpreal>(ctx, options, corrections);
}

RetCode ECHMET_CC createSolverContext(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	return createSolverContextInternal<ECHMETReal>(ctx, chemSystem);
}

RetCode ECHMET_CC createSolverContextHighPrecision(SolverContext *&ctx, const SysComp::ChemicalSystem &chemSystem) noexcept
{
	NumericPrecisionSetter<mpfr::mpreal> nps{};

	return createSolverContextInternal<mpfr::mpreal>(ctx, chemSystem);
}

InstructionSet detectInstructionSet() noexcept
{
	const CPUSIMD simd = cpuSupportedSIMD();

	if (simd.AVX512.F && simd.AVX512.DQ && simd.FMA3 && simd.AVX2 && simd.AVX && simd.SSE42 && simd.SSE41 && simd.SSSE3) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Using AVX512-optimized solver\n"));
		return InstructionSet::AVX512;
	}
	if (simd.FMA3 && simd.AVX2 && simd.AVX && simd.SSE42 && simd.SSE41 && simd.SSSE3) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Using FMA3-optimized solver\n"));
		return InstructionSet::FMA3;
	}
	if (simd.AVX && simd.SSE42 && simd.SSE41 && simd.SSSE3) {
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
