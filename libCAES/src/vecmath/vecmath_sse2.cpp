#include "vecmath.h"

#include <new>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif // ECHMET_COMPILER_

#define MK_VD2(v) { v, v }

namespace ECHMET {
namespace CAES {

typedef typename VecMath<InstructionSet::SSE2>::VD VD;

constexpr const uint64_t VecMath<InstructionSet::SSE2>::ZERO_BLOCK[VDType<InstructionSet::SSE2>::ALIGNMENT_BYTES / sizeof(double)] = MK_VD2(0);
constexpr const VD VecMath<InstructionSet::SSE2>::INFS = MK_VD2(VecMathCommon::DBL_INF);
constexpr const VD VecMath<InstructionSet::SSE2>::MINUS_INFS = MK_VD2(VecMathCommon::DBL_INF);
constexpr const VD VecMath<InstructionSet::SSE2>::MINUS_ONE = MK_VD2(-1.0);
constexpr const VD VecMath<InstructionSet::SSE2>::HALF = MK_VD2(0.5);
constexpr const VD VecMath<InstructionSet::SSE2>::ONE = MK_VD2(1.0);
constexpr const VD VecMath<InstructionSet::SSE2>::TWO = MK_VD2(2.0);
constexpr const VD VecMath<InstructionSet::SSE2>::MAXL10 = MK_VD2(VecMathCommon::MAXL10);
constexpr const VD VecMath<InstructionSet::SSE2>::MINUS_MAXL10 = MK_VD2(-VecMathCommon::MAXL10);

constexpr const VD VecMath<InstructionSet::SSE2>::LG102A = MK_VD2(VecMathCommon::LG102A);
constexpr const VD VecMath<InstructionSet::SSE2>::LG102B = MK_VD2(VecMathCommon::LG102B);
constexpr const VD VecMath<InstructionSet::SSE2>::LOG210 = MK_VD2(VecMathCommon::LOG210);

constexpr const VD VecMath<InstructionSet::SSE2>::L102A_LOG10 = MK_VD2(VecMathCommon::L102A_LOG10);
constexpr const VD VecMath<InstructionSet::SSE2>::L102B_LOG10 = MK_VD2(VecMathCommon::L102B_LOG10);
constexpr const VD VecMath<InstructionSet::SSE2>::L10EA_LOG10 = MK_VD2(VecMathCommon::L10EA_LOG10);
constexpr const VD VecMath<InstructionSet::SSE2>::L10EB_LOG10 = MK_VD2(VecMathCommon::L10EB_LOG10);
constexpr const VD VecMath<InstructionSet::SSE2>::SQRTH_LOG10 = MK_VD2(VecMathCommon::SQRTH_LOG10);

VecMath<InstructionSet::SSE2>::VecMath() :
	POne MK_VD2(*(double *)&VecMathCommon::PExp10[0]),
	PTwo MK_VD2(*(double *)&VecMathCommon::PExp10[4]),
	PThree MK_VD2(*(double *)&VecMathCommon::PExp10[8]),
	PFour MK_VD2(*(double *)&VecMathCommon::PExp10[12]),

	QOne MK_VD2(*(double *)&VecMathCommon::QExp10[0]),
	QTwo MK_VD2(*(double *)&VecMathCommon::QExp10[4]),
	QThree MK_VD2(*(double *)&VecMathCommon::QExp10[8]),

	POneLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[0]),
	PTwoLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[4]),
	PThreeLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[8]),
	PFourLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[12]),
	PFiveLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[16]),
	PSixLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[20]),
	PSevenLog10 MK_VD2(*(double*)&VecMathCommon::PLog10[24]),

	QOneLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[0]),
	QTwoLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[4]),
	QThreeLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[8]),
	QFourLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[12]),
	QFiveLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[16]),
	QSixLog10 MK_VD2(*(double*)&VecMathCommon::QLog10[20])
{}

typename VecMath<InstructionSet::SSE2>::TD VecMath<InstructionSet::SSE2>::exp10m(TD x) const
{
	x = _mm_mul_pd(x, M128D(MINUS_ONE));

	/* Over and underflow checks */
	uint64_t CHECKS[2] ECHMET_ALIGNED_16;
	__m128d tmp = _mm_cmpgt_pd(x, M128D(MAXL10));
	_mm_store_pd((double *)CHECKS, tmp);
	if (std::memcmp(CHECKS, ZERO_BLOCK, 2 * sizeof(uint64_t)))
		return M128D(INFS);

	tmp = _mm_cmplt_pd(x, M128D(MINUS_MAXL10));
	_mm_store_pd((double *)CHECKS, tmp);
	if (std::memcmp(CHECKS, ZERO_BLOCK, 2 * sizeof(uint64_t)))
		return M128D(MINUS_INFS);

	/* px = floor(LOG210 * x + 0.5) */
	__m128d px = M128D(LOG210);
	px = _mm_mul_pd(px, x);
	px = _mm_add_pd(px, M128D(HALF));
	/* Floor by conversion to int32 and back with truncation */
	__m64 mm0 = _mm_cvttpd_pi32(px);
	tmp = _mm_cvtpi32_pd(mm0);
	/* Safety net - subtract one if the result is greater than input */
	__m128d mask = _mm_cmpgt_pd(tmp, px);
	mask = _mm_and_pd(mask, M128D(ONE));
	px = _mm_sub_pd(tmp, mask);
	__m128d n = px;

	/* x = x - px * LG102A - px * LG102B */
	tmp = _mm_mul_pd(px, M128D(LG102A));
	x = _mm_sub_pd(x, tmp);
	tmp = _mm_mul_pd(px, M128D(LG102B));
	x = _mm_sub_pd(x, tmp);

	__m128d xx = _mm_mul_pd(x, x);
	/* polevl P */
	tmp = M128D(POne);
	tmp = _mm_mul_pd(tmp, xx);
	tmp = _mm_add_pd(tmp, M128D(PTwo));
	tmp = _mm_mul_pd(tmp, xx);
	tmp = _mm_add_pd(tmp, M128D(PThree));
	tmp = _mm_mul_pd(tmp, xx);
	tmp = _mm_add_pd(tmp, M128D(PFour));
	/* px = x * polevlP */
	px = _mm_mul_pd(x, tmp);

	/* polevl Q */
	tmp = M128D(QOne);
	tmp = _mm_mul_pd(tmp, xx);
	tmp = _mm_add_pd(tmp, M128D(QTwo));
	tmp = _mm_mul_pd(tmp, xx);
	tmp = _mm_add_pd(tmp, M128D(QThree));
	tmp = _mm_sub_pd(tmp, px); /* polevlQ - px */
	x = _mm_div_pd(px, tmp);
	x = _mm_mul_pd(x, M128D(TWO));
	x = _mm_add_pd(x, M128D(ONE));

	/*VD2 d_n;
	VD2 d_x;
	_mm_store_pd(d_n, n);
	_mm_store_pd(d_x, x);

	x = _mm_set_pd(cephes_ldexp(d_x[1], d_n[1]), cephes_ldexp(d_x[0], d_n[0]));*/
	x = _mm_set_pd(VecMathCommon::cephes_ldexp(x[1], n[1]), VecMathCommon::cephes_ldexp(x[0], n[0])); /* Array-like access to SIMD types is a GCC 4.6+ extenstion! */

	return x;
}

typename VecMath<InstructionSet::SSE2>::TD VecMath<InstructionSet::SSE2>::mlog10(const double *ECHMET_RESTRICT_PTR inx) const
{
	VD vx;
	int32_t ve[2] ECHMET_ALIGNED_16;
	vx[0] = VecMathCommon::cephes_frexp(inx[0], &ve[0]);
	vx[1] = VecMathCommon::cephes_frexp(inx[1], &ve[1]);

	__m128d x = M128D(vx);
	__m128d e = _mm_cvtpi32_pd(M64(ve));

	/* if (x < SQRTH)
	 * then
	 *	e = e - 1
	 *	x = 2x - 1
	 * else
	 *	x = x - 1
	 */
	__m128d tmp = _mm_cmplt_pd(x, M128D(SQRTH_LOG10));
	__m128d m = _mm_and_pd(tmp, M128D(ONE));
	e = _mm_sub_pd(e, m);	/* Conditional e = e - 1 */
	m = _mm_and_pd(tmp, x);
	x = _mm_add_pd(x, m);	/* Conditional x = x + x */
	x = _mm_sub_pd(x, M128D(ONE));

	__m128d xx = _mm_mul_pd(x, x);
	/* polevl P */
	tmp = M128D(POneLog10);
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PTwoLog10));
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PThreeLog10));
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PFourLog10));
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PFiveLog10));
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PSixLog10));
	tmp = _mm_mul_pd(tmp, x);
	tmp = _mm_add_pd(tmp, M128D(PSevenLog10));
	tmp = _mm_mul_pd(tmp, x);

	/* polevl Q */
	m = M128D(QOneLog10);
	m = _mm_mul_pd(m, x);
	m = _mm_add_pd(m, M128D(QTwoLog10));
	m = _mm_mul_pd(m, x);
	m = _mm_add_pd(m, M128D(QThreeLog10));
	m = _mm_mul_pd(m, x);
	m = _mm_add_pd(m, M128D(QFourLog10));
	m = _mm_mul_pd(m, x);
	m = _mm_add_pd(m, M128D(QFiveLog10));
	m = _mm_mul_pd(m, x);
	m = _mm_add_pd(m, M128D(QSixLog10));

	m = _mm_div_pd(tmp, m);
	m = _mm_mul_pd(xx, m);
	tmp = _mm_mul_pd(xx, M128D(HALF));
	m = _mm_sub_pd(m, tmp);

	/* xx = (x + y) * L10EB */
	tmp = _mm_add_pd(x, m);
	xx = _mm_mul_pd(tmp, M128D(L10EB_LOG10));

	/* xx = xx + y * L10EA */
	tmp = _mm_mul_pd(m, M128D(L10EA_LOG10));
	xx = _mm_add_pd(xx, tmp);

	/* xx = xx + x * L10EA */
	tmp = _mm_mul_pd(x, M128D(L10EA_LOG10));
	xx = _mm_add_pd(xx, tmp);

	/* xx = xx + e * L102B */
	tmp = _mm_mul_pd(e, M128D(L102B_LOG10));
	xx = _mm_add_pd(xx, tmp);

	/* xx = xx + e * L102A */
	tmp = _mm_mul_pd(e, M128D(L102A_LOG10));
	xx = _mm_add_pd(xx, tmp);

	xx = _mm_mul_pd(xx, M128D(MINUS_ONE));

	return xx;
}

double VecMath<InstructionSet::SSE2>::exp10m_single(const double x) noexcept
{
	return VecMathCommon::exp10m_single(x);
}

double VecMath<InstructionSet::SSE2>::mlog10_single(const double x) noexcept
{
	return VecMathCommon::mlog10_single(x);
}

void VecMath<InstructionSet::SSE2>::operator delete(void *ptr)
{
	_mm_free(ptr);
}

void * VecMath<InstructionSet::SSE2>::operator new(const size_t sz)
{
	void *ptr = _mm_malloc(sz, VDType<InstructionSet::SSE2>::ALIGNMENT_BYTES);
	if (ptr == nullptr)
		throw std::bad_alloc{};

	return ptr;
}

} // namespace CAES
} // namespace ECHMET
