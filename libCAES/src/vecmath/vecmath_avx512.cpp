#include "vecmath.h"

#include <new>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif // ECHMET_COMPILER_

#define MK_VD8(v) { v, v, v, v, v, v, v, v }

namespace ECHMET {
namespace CAES {

typedef typename VecMath<InstructionSet::AVX512>::VD VD;

const VD VecMath<InstructionSet::AVX512>::INFS = MK_VD8(VecMathCommon::DBL_INF);
const VD VecMath<InstructionSet::AVX512>::MINUS_INFS = MK_VD8(VecMathCommon::DBL_INF);
const VD VecMath<InstructionSet::AVX512>::MINUS_ONE = MK_VD8(-1.0);
const VD VecMath<InstructionSet::AVX512>::HALF = MK_VD8(0.5);
const VD VecMath<InstructionSet::AVX512>::ONE = MK_VD8(1.0);
const VD VecMath<InstructionSet::AVX512>::TWO = MK_VD8(2.0);
const VD VecMath<InstructionSet::AVX512>::MAXL10 = MK_VD8(VecMathCommon::MAXL10);
const VD VecMath<InstructionSet::AVX512>::MINUS_MAXL10 = MK_VD8(-VecMathCommon::MAXL10);
const VD VecMath<InstructionSet::AVX512>::LG102A = MK_VD8(VecMathCommon::LG102A);
const VD VecMath<InstructionSet::AVX512>::LG102B = MK_VD8(VecMathCommon::LG102B);
const VD VecMath<InstructionSet::AVX512>::LOG210 = MK_VD8(VecMathCommon::LOG210);
const VD VecMath<InstructionSet::AVX512>::L102A_LOG10 = MK_VD8(VecMathCommon::L102A_LOG10);
const VD VecMath<InstructionSet::AVX512>::L102B_LOG10 = MK_VD8(VecMathCommon::L102B_LOG10);
const VD VecMath<InstructionSet::AVX512>::L10EA_LOG10 = MK_VD8(VecMathCommon::L10EA_LOG10);
const VD VecMath<InstructionSet::AVX512>::L10EB_LOG10 = MK_VD8(VecMathCommon::L10EB_LOG10);
const VD VecMath<InstructionSet::AVX512>::SQRTH_LOG10 = MK_VD8(VecMathCommon::SQRTH_LOG10);

VecMath<InstructionSet::AVX512>::VecMath() :
	POne MK_VD8(*(double *)&VecMathCommon::PExp10[0]),
	PTwo MK_VD8(*(double *)&VecMathCommon::PExp10[4]),
	PThree MK_VD8(*(double *)&VecMathCommon::PExp10[8]),
	PFour MK_VD8(*(double *)&VecMathCommon::PExp10[12]),

	QOne MK_VD8(*(double *)&VecMathCommon::QExp10[0]),
	QTwo MK_VD8(*(double *)&VecMathCommon::QExp10[4]),
	QThree MK_VD8(*(double *)&VecMathCommon::QExp10[8]),

	POneLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[0]),
	PTwoLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[4]),
	PThreeLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[8]),
	PFourLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[12]),
	PFiveLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[16]),
	PSixLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[20]),
	PSevenLog10 MK_VD8(*(double*)&VecMathCommon::PLog10[24]),

	QOneLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[0]),
	QTwoLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[4]),
	QThreeLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[8]),
	QFourLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[12]),
	QFiveLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[16]),
	QSixLog10 MK_VD8(*(double*)&VecMathCommon::QLog10[20])
{}

#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
void VecMath<InstructionSet::AVX512>::exp10m(const double *ECHMET_RESTRICT_PTR inx, double *ECHMET_RESTRICT_PTR outx) const
#else
typename VecMath<InstructionSet::AVX512>::TD VecMath<InstructionSet::AVX512>::exp10m(TD x) const
#endif
{
#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
	__m512d x = _mm512_loadu_pd(inx);
#endif
	x = _mm512_mul_pd(x, M512D(MINUS_ONE));

	/* Over and underflow checks */
	 __mmask8 mask = _mm512_cmp_pd_mask(x, M512D(MAXL10), _CMP_GT_OS);
	if (mask) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
		_mm512_store_pd(outx, M512D(INFS));
		return;
	#else
		return M512D(INFS);
	#endif
	}

	mask = _mm512_cmp_pd_mask(x, M512D(MINUS_MAXL10), _CMP_LT_OS);
	if (mask) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
		_mm512_store_pd(outx, M512D(MINUS_INFS));
		return;
	#else
		return M512D(MINUS_INFS);
	#endif
	}

	/* px = floor(LOG210 * x + 0.5) */
	__m512d px = M512D(LOG210);
	px = _mm512_mul_pd(px, x);
	px = _mm512_add_pd(px, M512D(HALF));
	/* Floor by conversion to int32 and back with truncation */
	__m256i emm0 = _mm512_cvttpd_epi32(px);
	__m512d tmp = _mm512_cvtepi32_pd(emm0);
	/* Safety net - subtract one if the result is greater than input */
	mask = _mm512_cmp_pd_mask(tmp, px, _CMP_GT_OS);
	px = _mm512_mask_sub_pd(tmp, mask, tmp, M512D(ONE));
	__m512d n = px;

	/* x = x - px * LG102A - px * LG102B */
	tmp = _mm512_mul_pd(px, M512D(LG102A));
	x = _mm512_sub_pd(x, tmp);
	tmp = _mm512_mul_pd(px, M512D(LG102B));
	x = _mm512_sub_pd(x, tmp);

	__m512d xx = _mm512_mul_pd(x, x);
	/* polevl P */
	tmp = M512D(POne);
	tmp = _mm512_fmadd_pd(tmp, xx, M512D(PTwo));
	tmp = _mm512_fmadd_pd(tmp, xx, M512D(PThree));
	tmp = _mm512_fmadd_pd(tmp, xx, M512D(PFour));
	/* px = x * polevlP */
	px = _mm512_mul_pd(x, tmp);

	/* polevl Q */
	tmp = M512D(QOne);
	tmp = _mm512_fmadd_pd(tmp, xx, M512D(QTwo));
	tmp = _mm512_fmadd_pd(tmp, xx, M512D(QThree));
	tmp = _mm512_sub_pd(tmp, px); /* polevlQ - px */
	x = _mm512_div_pd(px, tmp);
	x = _mm512_fmadd_pd(x, M512D(TWO), M512D(ONE));


#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	x = _mm512_set_pd(VecMathCommon::cephes_ldexp(x[7], n[7]),
			  VecMathCommon::cephes_ldexp(x[6], n[6]),
			  VecMathCommon::cephes_ldexp(x[5], n[5]),
			  VecMathCommon::cephes_ldexp(x[4], n[4]),
			  VecMathCommon::cephes_ldexp(x[3], n[3]),
			  VecMathCommon::cephes_ldexp(x[2], n[2]),
			  VecMathCommon::cephes_ldexp(x[1], n[1]),
			  VecMathCommon::cephes_ldexp(x[0], n[0])); /* Array-like access to SIMD types is a GCC 4.6+ extenstion! */
#else
	VD d_n;
	VD d_x;
	_mm512_store_pd(d_n, n);
	_mm512_store_pd(d_x, x);

	x = _mm512_set_pd(VecMathCommon::cephes_ldexp(d_x[7], d_n[7]),
			  VecMathCommon::cephes_ldexp(d_x[6], d_n[6]),
			  VecMathCommon::cephes_ldexp(d_x[5], d_n[5]),
			  VecMathCommon::cephes_ldexp(d_x[4], d_n[4]),
			  VecMathCommon::cephes_ldexp(d_x[3], d_n[3]),
			  VecMathCommon::cephes_ldexp(d_x[2], d_n[2]),
			  VecMathCommon::cephes_ldexp(d_x[1], d_n[1]),
			  VecMathCommon::cephes_ldexp(d_x[0], d_n[0]));
#endif // ECHMET_COMPILER_


#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
	_mm512_storeu_pd(outx, x);
#else
	return x;
#endif
}

typename VecMath<InstructionSet::AVX512>::TD VecMath<InstructionSet::AVX512>::mlog10(const double *ECHMET_RESTRICT_PTR inx) const
{
	VD vx;
	int32_t ECHMET_ALIGNED_BEF_64 ve[8] ECHMET_ALIGNED_AFT_64;
	vx[0] = VecMathCommon::cephes_frexp(inx[0], &ve[0]);
	vx[1] = VecMathCommon::cephes_frexp(inx[1], &ve[1]);
	vx[2] = VecMathCommon::cephes_frexp(inx[2], &ve[2]);
	vx[3] = VecMathCommon::cephes_frexp(inx[3], &ve[3]);
	vx[4] = VecMathCommon::cephes_frexp(inx[4], &ve[4]);
	vx[5] = VecMathCommon::cephes_frexp(inx[5], &ve[5]);
	vx[6] = VecMathCommon::cephes_frexp(inx[6], &ve[6]);
	vx[7] = VecMathCommon::cephes_frexp(inx[7], &ve[7]);

	__m512d x = M512D(vx);
	__m512d e = _mm512_cvtepi32_pd(M256I(ve));

	/* if (x < SQRTH)
	 * then
	 *	e = e - 1
	 *	x = 2x - 1
	 * else
	 *	x = x - 1
	 */
	__mmask8 mask = _mm512_cmp_pd_mask(x, M512D(SQRTH_LOG10), _CMP_LT_OS);
	e = _mm512_mask_sub_pd(e, mask, e, M512D(ONE));	/* Conditional e = e - 1 */
	x = _mm512_mask_add_pd(x, mask, x, M512D(ONE));	/* Conditional x = x + x */
	x = _mm512_sub_pd(x, M512D(ONE));

	__m512d xx = _mm512_mul_pd(x, x);
	/* polevl P */
	__m512d tmp = M512D(POneLog10);
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PTwoLog10));
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PThreeLog10));
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PFourLog10));
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PFiveLog10));
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PSixLog10));
	tmp = _mm512_fmadd_pd(tmp, x, M512D(PSevenLog10));
	tmp = _mm512_mul_pd(tmp, x);

	/* polevl Q */
	__m512d m = M512D(QOneLog10);
	m = _mm512_fmadd_pd(m, x, M512D(QTwoLog10));
	m = _mm512_fmadd_pd(m, x, M512D(QThreeLog10));
	m = _mm512_fmadd_pd(m, x, M512D(QFourLog10));
	m = _mm512_fmadd_pd(m, x, M512D(QFiveLog10));
	m = _mm512_fmadd_pd(m, x, M512D(QSixLog10));

	m = _mm512_div_pd(tmp, m);
	m = _mm512_mul_pd(xx, m);
	tmp = _mm512_mul_pd(xx, M512D(HALF));
	m = _mm512_sub_pd(m, tmp);

	/* xx = (x + y) * L10EB */
	tmp = _mm512_add_pd(x, m);
	xx = _mm512_mul_pd(tmp, M512D(L10EB_LOG10));

	/* xx = xx + y * L10EA */
	xx = _mm512_fmadd_pd(m, M512D(L10EA_LOG10), xx);

	/* xx = xx + x * L10EA */
	xx = _mm512_fmadd_pd(x, M512D(L10EA_LOG10), xx);

	/* xx = xx + e * L102B */
	xx = _mm512_fmadd_pd(e, M512D(L102B_LOG10), xx);

	/* xx = xx + e * L102A */
	xx = _mm512_fmadd_pd(e, M512D(L102A_LOG10), xx);

	xx = _mm512_mul_pd(xx, M512D(MINUS_ONE));

	return xx;
}

double VecMath<InstructionSet::AVX512>::exp10m_single(const double x) noexcept
{
	return VecMathCommon::exp10m_single(x);
}

double VecMath<InstructionSet::AVX512>::mlog10_single(const double x) noexcept
{
	return VecMathCommon::mlog10_single(x);
}

void VecMath<InstructionSet::AVX512>::operator delete(void *ptr)
{
	_mm_free(ptr);
}

void * VecMath<InstructionSet::AVX512>::operator new(const size_t sz)
{
	void *ptr = _mm_malloc(sz, VDType<InstructionSet::AVX512>::ALIGNMENT_BYTES);
	if (ptr == nullptr)
		throw std::bad_alloc{};

	return ptr;
}

} // namespace CAES
} // namespace ECHMET
