#include "vecmath/vecmath.h"

#define ECHMET_IMPORT_INTERNAL
#include <echmetelems.h>

#ifdef ECHMET_DEBUG_CODE
#include <cstdio>
#endif // ECHMET_DEBUG_CODE

namespace ECHMET {
namespace CAES {

InstructionSet detectInstructionSet()
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

} // namespace CAES
} // namespace ECHMET
