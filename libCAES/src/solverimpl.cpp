#include "vecmath/vecmath.h"
#include "cpufeatures.h"

#include <echmetelems.h>

namespace ECHMET {
namespace CAES {

InstructionSet detectInstructionSet()
{
	const auto simd = CPUFeatures::SIMD();

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

	return InstructionSet::GENERIC;
}

} // namespace CAES
} // namespace ECHMET
