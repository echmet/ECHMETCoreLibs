#ifndef ECHMET_CAES_VECMATH_H
#define ECHMET_CAES_VECMATH_H

#include <cstring>
#include <cstdint>
#include <limits>
#include <new>
#include "../veccommon.h"

namespace ECHMET {
namespace CAES {

class VecMathCommon {
public:
	/* Constants to compute 10^x borrowed from Cephes math library */
	constexpr static const double LOG210 = 3.32192809488736234787e0;
	constexpr static const double LG102A = 3.01025390625000000000e-1;
	constexpr static const double LG102B = 4.60503898119521373889e-6;
	constexpr static const double DBL_INF = std::numeric_limits<double>::infinity();
	constexpr static const double MAXL10 = 308.2547155599167;

	/* Constants to compute log10 borrowed from Cephes math library */
	constexpr static const double L102A_LOG10 = 3.0078125E-1;
	constexpr static const double L102B_LOG10 = 2.48745663981195213739E-4;
	constexpr static const double L10EA_LOG10 = 4.3359375E-1;
	constexpr static const double L10EB_LOG10 = 7.00731903251827651129E-4;
	constexpr static const double SQRTH_LOG10 = 0.70710678118654752440;

	static const unsigned short PExp10[];
	static const unsigned short QExp10[];
	static const unsigned short PLog10[];
	static const unsigned short QLog10[];

	static double cephes_frexp(const double x, int32_t *pw2) noexcept;
	static double cephes_ldexp(const double x, int32_t pw2) noexcept;

	static double exp10m_single(double x) noexcept;
	static double mlog10_single(double x) noexcept;
};

template <>
class VDType<InstructionSet::SSE2> {
public:
	typedef __m128d TD;
	typedef double ECHMET_ALIGNED_BEF_16 VD[2] ECHMET_ALIGNED_AFT_16;

	constexpr static const size_t ALIGNMENT_BYTES = 16;
};

template <>
class VDType<InstructionSet::AVX> {
public:
	typedef __m256d TD;
	typedef double ECHMET_ALIGNED_BEF_32 VD[4] ECHMET_ALIGNED_AFT_32;

	constexpr static const size_t ALIGNMENT_BYTES = 32;
};

template <>
class VDType<InstructionSet::FMA3> {
public:
	typedef __m256d TD;
	typedef double ECHMET_ALIGNED_BEF_32 VD[4] ECHMET_ALIGNED_AFT_32;

	constexpr static const size_t ALIGNMENT_BYTES = 32;
};

#ifndef ECHMET_DISABLE_AVX512

template <>
class VDType<InstructionSet::AVX512> {
public:
	typedef __m512d TD;
	typedef double ECHMET_ALIGNED_BEF_64 VD[8] ECHMET_ALIGNED_AFT_64;

	constexpr static const size_t ALIGNMENT_BYTES = 64;
};

#endif // ECHMET_DISABLE_AVX512

#define MAKE_VECMATH(iset) \
	template <> \
	class VecMath<InstructionSet::iset> { \
	public: \
		typedef typename VDType<InstructionSet::iset>::TD TD; \
		typedef typename VDType<InstructionSet::iset>::VD VD; \
 \
		VecMath(); \
		TD exp10m(TD x) const; \
		static double exp10m_single(double x) noexcept; \
		TD mlog10(const double *ECHMET_RESTRICT_PTR inx) const; \
		static double mlog10_single(double x) noexcept; \
 \
		void operator delete(void *ptr); \
		void * operator new(const size_t sz); \
 \
		static const uint64_t ZERO_BLOCK[]; \
		static const VD INFS; \
		static const VD MINUS_INFS; \
		static const VD MINUS_ONE; \
		static const VD HALF; \
		static const VD ONE; \
		static const VD TWO; \
		static const VD MAXL10; \
		static const VD MINUS_MAXL10; \
 \
		static const VD LG102A; \
		static const VD LG102B; \
		static const VD LOG210; \
 \
		static const VD L102A_LOG10; \
		static const VD L102B_LOG10; \
		static const VD L10EA_LOG10; \
		static const VD L10EB_LOG10; \
		static const VD SQRTH_LOG10; \
 \
		VD POne; \
		VD PTwo; \
		VD PThree; \
		VD PFour; \
		VD QOne; \
		VD QTwo; \
		VD QThree; \
 \
		VD POneLog10; \
		VD PTwoLog10; \
		VD PThreeLog10; \
		VD PFourLog10; \
		VD PFiveLog10; \
		VD PSixLog10; \
		VD PSevenLog10; \
 \
		VD QOneLog10; \
		VD QTwoLog10; \
		VD QThreeLog10; \
		VD QFourLog10; \
		VD QFiveLog10; \
		VD QSixLog10; \
	}

MAKE_VECMATH(SSE2);
MAKE_VECMATH(AVX);
MAKE_VECMATH(FMA3);
#ifndef ECHMET_DISABLE_AVX512
MAKE_VECMATH(AVX512);
#endif // ECHMET_DISABLE_AVX512

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_VECMATH_H
