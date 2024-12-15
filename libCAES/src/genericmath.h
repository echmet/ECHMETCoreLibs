#ifndef ECHMET_CAES_GENERICMATH_H
#define ECHMET_CAES_GENERICMATH_H

#include <cstdint>

namespace ECHMET {
namespace CAES {

enum InstructionSet {
	GENERIC,
	SSE2,
	AVX,
	FMA3,
	AVX512
};

template <InstructionSet ISet>
class VDType {
public:
	constexpr static const size_t ALIGNMENT_BYTES = 16;
};

template <InstructionSet ISet>
class VecMath {
};


/*
 * TODO: There should be a generic ifdef flag to enable this allocator
 * if there are no architecture-specific allocators available.
 * We currently implement architecture-specific stuff only for x86
 * so using the _X86 flag below is okay
 */
#ifndef ECHMET_USE_X86_EXTENSIONS

template <typename T, size_t Alignment>
class AlignedAllocator {
public:
	static T * alloc(const size_t N)
	{
		return new T[N];
	}

	static void free(T *ptr)
	{
		return delete[] ptr;
	}
};

#endif // ECHMET_USE_X86_EXTENSIONS

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_GENERICMATH_H
