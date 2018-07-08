#include "solverinternal.h"
#include "vecmath/vecmath.h"

namespace ECHMET {
namespace CAES {

template <> template <>
void SolverInternal<double, InstructionSet::SSE2>::VectorizedDelogifier<InstructionSet::SSE2>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "SSE2 delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m128d _s = M128D(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}

template <> template <>
void SolverInternal<double, InstructionSet::SSE2>::VectorizedLogifier<InstructionSet::SSE2>::operator()(double  *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "SSE2 logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m128d _s = m_vecMath.mlog10(src + idx);

		_mm_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}

} // namespace CAES
} // namespace ECHMET
