#include "vecmath.h"

#include <new>

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif // ECHMET_COMPILER_

#define MK_VD4(v) { v, v, v, v }

namespace ECHMET {
namespace CAES {

typedef typename VecMath<InstructionSet::FMA3>::VD VD;

constexpr const uint64_t VecMath<InstructionSet::FMA3>::ZERO_BLOCK[VDType<InstructionSet::FMA3>::ALIGNMENT_BYTES / sizeof(double)] = MK_VD4(0);
constexpr const VD VecMath<InstructionSet::FMA3>::INFS = MK_VD4(VecMathCommon::DBL_INF);
constexpr const VD VecMath<InstructionSet::FMA3>::MINUS_INFS = MK_VD4(VecMathCommon::DBL_INF);
constexpr const VD VecMath<InstructionSet::FMA3>::MINUS_ONE = MK_VD4(-1.0);
constexpr const VD VecMath<InstructionSet::FMA3>::HALF = MK_VD4(0.5);
constexpr const VD VecMath<InstructionSet::FMA3>::ONE = MK_VD4(1.0);
constexpr const VD VecMath<InstructionSet::FMA3>::TWO = MK_VD4(2.0);
constexpr const VD VecMath<InstructionSet::FMA3>::MAXL10 = MK_VD4(VecMathCommon::MAXL10);
constexpr const VD VecMath<InstructionSet::FMA3>::MINUS_MAXL10 = MK_VD4(-VecMathCommon::MAXL10);
constexpr const VD VecMath<InstructionSet::FMA3>::LG102A = MK_VD4(VecMathCommon::LG102A);
constexpr const VD VecMath<InstructionSet::FMA3>::LG102B = MK_VD4(VecMathCommon::LG102B);
constexpr const VD VecMath<InstructionSet::FMA3>::LOG210 = MK_VD4(VecMathCommon::LOG210);
constexpr const VD VecMath<InstructionSet::FMA3>::L102A_LOG10 = MK_VD4(VecMathCommon::L102A_LOG10);
constexpr const VD VecMath<InstructionSet::FMA3>::L102B_LOG10 = MK_VD4(VecMathCommon::L102B_LOG10);
constexpr const VD VecMath<InstructionSet::FMA3>::L10EA_LOG10 = MK_VD4(VecMathCommon::L10EA_LOG10);
constexpr const VD VecMath<InstructionSet::FMA3>::L10EB_LOG10 = MK_VD4(VecMathCommon::L10EB_LOG10);
constexpr const VD VecMath<InstructionSet::FMA3>::SQRTH_LOG10 = MK_VD4(VecMathCommon::SQRTH_LOG10);

VecMath<InstructionSet::FMA3>::VecMath() :
	POne MK_VD4(*(double *)&VecMathCommon::PExp10[0]),
	PTwo MK_VD4(*(double *)&VecMathCommon::PExp10[4]),
	PThree MK_VD4(*(double *)&VecMathCommon::PExp10[8]),
	PFour MK_VD4(*(double *)&VecMathCommon::PExp10[12]),

	QOne MK_VD4(*(double *)&VecMathCommon::QExp10[0]),
	QTwo MK_VD4(*(double *)&VecMathCommon::QExp10[4]),
	QThree MK_VD4(*(double *)&VecMathCommon::QExp10[8]),

	POneLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[0]),
	PTwoLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[4]),
	PThreeLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[8]),
	PFourLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[12]),
	PFiveLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[16]),
	PSixLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[20]),
	PSevenLog10 MK_VD4(*(double*)&VecMathCommon::PLog10[24]),

	QOneLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[0]),
	QTwoLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[4]),
	QThreeLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[8]),
	QFourLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[12]),
	QFiveLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[16]),
	QSixLog10 MK_VD4(*(double*)&VecMathCommon::QLog10[20])
{}

#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
void VecMath<InstructionSet::FMA3>::exp10m(const double *ECHMET_RESTRICT_PTR inx, double *ECHMET_RESTRICT_PTR outx) const
#else
typename VecMath<InstructionSet::FMA3>::TD VecMath<InstructionSet::FMA3>::exp10m(TD x) const
#endif
{
#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
	__m256d x = _mm256_loadu_pd(inx);
#endif
	x = _mm256_mul_pd(x, M256D(MINUS_ONE));

	/* Over and underflow checks */
	uint64_t CHECKS[4] ECHMET_ALIGNED_32;
	__m256d tmp = _mm256_cmp_pd(x, M256D(MAXL10), _CMP_GT_OS);
	_mm256_store_pd((double *)CHECKS, tmp);
	if (std::memcmp(CHECKS, ZERO_BLOCK, 4 * sizeof(uint64_t))) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
		_mm256_store_pd(outx, M256D(INFS));
		return;
	#else
		return M256D(INFS);
	#endif
	}

	tmp = _mm256_cmp_pd(x, M256D(MINUS_MAXL10), _CMP_LT_OS);
	_mm256_store_pd((double *)CHECKS, tmp);
	if (std::memcmp(CHECKS, ZERO_BLOCK, 4 * sizeof(uint64_t))) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
		_mm256_store_pd(outx, M256D(MINUS_INFS));
		return;
	#else
		return M256D(MINUS_INFS);
	#endif
	}

	/* px = floor(LOG210 * x + 0.5) */
	__m256d px = M256D(LOG210);
	px = _mm256_mul_pd(px, x);
	px = _mm256_add_pd(px, M256D(HALF));
	/* Floor by conversion to int32 and back with truncation */
	__m128i emm0 = _mm256_cvttpd_epi32(px);
	tmp = _mm256_cvtepi32_pd(emm0);
	/* Safety net - subtract one if the result is greater than input */
	__m256d mask = _mm256_cmp_pd(tmp, px, _CMP_GT_OS);
	mask = _mm256_and_pd(mask, M256D(ONE));
	px = _mm256_sub_pd(tmp, mask);
	__m256d n = px;

	/* x = x - px * LG102A - px * LG102B */
	tmp = _mm256_mul_pd(px, M256D(LG102A));
	x = _mm256_sub_pd(x, tmp);
	tmp = _mm256_mul_pd(px, M256D(LG102B));
	x = _mm256_sub_pd(x, tmp);

	__m256d xx = _mm256_mul_pd(x, x);
	/* polevl P */
	tmp = M256D(POne);
	tmp = _mm256_fmadd_pd(tmp, xx, M256D(PTwo));
	tmp = _mm256_fmadd_pd(tmp, xx, M256D(PThree));
	tmp = _mm256_fmadd_pd(tmp, xx, M256D(PFour));
	/* px = x * polevlP */
	px = _mm256_mul_pd(x, tmp);

	/* polevl Q */
	tmp = M256D(QOne);
	tmp = _mm256_fmadd_pd(tmp, xx, M256D(QTwo));
	tmp = _mm256_fmadd_pd(tmp, xx, M256D(QThree));
	tmp = _mm256_sub_pd(tmp, px); /* polevlQ - px */
	x = _mm256_div_pd(px, tmp);
	x = _mm256_fmadd_pd(x, M256D(TWO), M256D(ONE));

	/*VD2 d_n;
	VD2 d_x;
	_mm_store_pd(d_n, n);
	_mm_store_pd(d_x, x);

	x = _mm_set_pd(cephes_ldexp(d_x[1], d_n[1]), cephes_ldexp(d_x[0], d_n[0]));*/
	x = _mm256_set_pd(VecMathCommon::cephes_ldexp(x[3], n[3]),
			  VecMathCommon::cephes_ldexp(x[2], n[2]),
			  VecMathCommon::cephes_ldexp(x[1], n[1]),
			  VecMathCommon::cephes_ldexp(x[0], n[0])); /* Array-like access to SIMD types is a GCC 4.6+ extenstion! */

#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__)
	_mm256_storeu_pd(outx, x);
#else
	return x;
#endif
}

typename VecMath<InstructionSet::FMA3>::TD VecMath<InstructionSet::FMA3>::mlog10(const double *ECHMET_RESTRICT_PTR inx) const
{
	VD vx;
	int32_t ve[4] ECHMET_ALIGNED_32;
	vx[0] = VecMathCommon::cephes_frexp(inx[0], &ve[0]);
	vx[1] = VecMathCommon::cephes_frexp(inx[1], &ve[1]);
	vx[2] = VecMathCommon::cephes_frexp(inx[2], &ve[2]);
	vx[3] = VecMathCommon::cephes_frexp(inx[3], &ve[3]);

	__m256d x = M256D(vx);
	__m256d e = _mm256_cvtepi32_pd(M128I(ve));

	/* if (x < SQRTH)
	 * then
	 *	e = e - 1
	 *	x = 2x - 1
	 * else
	 *	x = x - 1
	 */
	__m256d tmp = _mm256_cmp_pd(x, M256D(SQRTH_LOG10), _CMP_LT_OS);
	__m256d m = _mm256_and_pd(tmp, M256D(ONE));
	e = _mm256_sub_pd(e, m);	/* Conditional e = e - 1 */
	m = _mm256_and_pd(tmp, x);
	x = _mm256_add_pd(x, m);	/* Conditional x = x + x */
	x = _mm256_sub_pd(x, M256D(ONE));

	__m256d xx = _mm256_mul_pd(x, x);
	/* polevl P */
	tmp = M256D(POneLog10);
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PTwoLog10));
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PThreeLog10));
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PFourLog10));
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PFiveLog10));
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PSixLog10));
	tmp = _mm256_fmadd_pd(tmp, x, M256D(PSevenLog10));
	tmp = _mm256_mul_pd(tmp, x);

	/* polevl Q */
	m = M256D(QOneLog10);
	m = _mm256_fmadd_pd(m, x, M256D(QTwoLog10));
	m = _mm256_fmadd_pd(m, x, M256D(QThreeLog10));
	m = _mm256_fmadd_pd(m, x, M256D(QFourLog10));
	m = _mm256_fmadd_pd(m, x, M256D(QFiveLog10));
	m = _mm256_fmadd_pd(m, x, M256D(QSixLog10));

	m = _mm256_div_pd(tmp, m);
	m = _mm256_mul_pd(xx, m);
	tmp = _mm256_mul_pd(xx, M256D(HALF));
	m = _mm256_sub_pd(m, tmp);

	/* xx = (x + y) * L10EB */
	tmp = _mm256_add_pd(x, m);
	xx = _mm256_mul_pd(tmp, M256D(L10EB_LOG10));

	/* xx = xx + y * L10EA */
	xx = _mm256_fmadd_pd(m, M256D(L10EA_LOG10), xx);

	/* xx = xx + x * L10EA */
	xx = _mm256_fmadd_pd(x, M256D(L10EA_LOG10), xx);

	/* xx = xx + e * L102B */
	xx = _mm256_fmadd_pd(e, M256D(L102B_LOG10), xx);

	/* xx = xx + e * L102A */
	xx = _mm256_fmadd_pd(e, M256D(L102A_LOG10), xx);

	xx = _mm256_mul_pd(xx, M256D(MINUS_ONE));

	return xx;
}

double VecMath<InstructionSet::FMA3>::exp10m_single(const double x) noexcept
{
	return VecMathCommon::exp10m_single(x);
}

double VecMath<InstructionSet::FMA3>::mlog10_single(const double x) noexcept
{
	return VecMathCommon::mlog10_single(x);
}

void VecMath<InstructionSet::FMA3>::operator delete(void *ptr)
{
	_mm_free(ptr);
}

void * VecMath<InstructionSet::FMA3>::operator new(const size_t sz)
{
	void *ptr = _mm_malloc(sz, VDType<InstructionSet::FMA3>::ALIGNMENT_BYTES);
	if (ptr == nullptr)
		throw std::bad_alloc{};

	return ptr;
}

} // namespace CAES
} // namespace ECHMET
