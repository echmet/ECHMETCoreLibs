#ifndef ECHMET_CAES_VECCOMMON_H
#define ECHMET_CAES_VECCOMMON_H

#include <new>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	#include <x86intrin.h>
#else
	#include <xmmintrin.h>
	#include <immintrin.h>
	#include <emmintrin.h>
#endif // ECHMET_COMPILER_

#define ECHMET_IMPORT_INTERNAL
	#include <echmetmodule.h>
#undef ECHMET_IMPORT_INTERNAL

#define M64(v) *(__m64 *)(v)
#define M128D(v) *(__m128d *)(v)
#define M128I(v) *(__m128i *)(v)
#define M256D(v) *(__m256d *)(v)

namespace ECHMET {
namespace CAES {

enum InstructionSet {
	GENERIC,
	SSE2,
	AVX,
	FMA3
};

template <typename T, size_t Alignment, bool RawAllocation>
class AlignedAllocatorWorker;

template <typename T, size_t Alignment>
class AlignedAllocatorWorker<T, Alignment, true> {
public:
	static T * alloc(const size_t N)
	{
		T *p = static_cast<T *>(_mm_malloc(sizeof(T) * N, Alignment));
		if (p == nullptr)
			throw std::bad_alloc{};

		return p;
	}

	static void free(T *ptr)
	{
		_mm_free(ptr);
	}
};

template <typename T, size_t Alignment>
class AlignedAllocatorWorker<T, Alignment, false> {
public:
	static T * alloc(const size_t N)
	{
		return new T[N];
	}

	static void free(T *ptr)
	{
		delete [] ptr;
	}
};

template <typename T, size_t Alignment>
class AlignedAllocator {
public:
	static T * alloc(const size_t N)
	{
		return AlignedAllocatorWorker<T, Alignment, std::is_fundamental<T>::value>::alloc(N);
	}

	static void free(T *ptr)
	{
		return AlignedAllocatorWorker<T, Alignment, std::is_fundamental<T>::value>::free(ptr);
	}
};

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_VECCOMMON_H
