#include "solverinternal.h"
#include "vecmath/vecmath.h"

namespace ECHMET {
namespace CAES {

template <> template <>
void SolverInternal<double, InstructionSet::AVX>::VectorizedDelogifier<InstructionSet::AVX>::operator()(double  *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
		m_vecMath.exp10m(src + idx, dst + idx);
	#else
		__m256d _s = M256D(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm256_store_pd(dst + idx, _s);
	#endif
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}

template <> template <>
void SolverInternal<double, InstructionSet::AVX>::VectorizedLogifier<InstructionSet::AVX>::operator()(double  *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m256d _s = m_vecMath.mlog10(src + idx);

		_mm256_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}

} // namespace CAES
} // namespace ECHMET
