#include "../containers/echmetfixedstring_p.h"
#include "../containers/echmetvec_p.h"
#include "cpufeatures.h"
#include "echmetelems.h"

#include <cstring>
#include <typeinfo>

#define _STRINGIFY(input) #input
#define ERROR_CODE_CASE(erCase) case RetCode::erCase: return _STRINGIFY(erCase)

#define SET_SIMD_FLAG(pub, intr, flag) pub.flag = intr.flag

namespace ECHMET {

RealVec * ECHMET_CC createRealVec(const size_t reserve) noexcept
{
	return createECHMETVec<ECHMETReal, false>(reserve);
}

RealVec * ECHMET_CC createRealVecPreallocated(const size_t size, ECHMETReal *const pool)
{
	return createECHMETVecPreallocated<ECHMETReal, false>(size, pool);
}

FixedString * ECHMET_CC createFixedString(const char *str) noexcept
{
	return new (std::nothrow) FixedStringImpl(str);
}

double ECHMET_CC ECHMETRealToDouble(const ECHMETReal &real) noexcept
{
#ifdef ECHMET_USE_HIGH_PRECISION
	return real.toDouble();
#else
	return real;
#endif // ECHMET_USE_HIGH_PRECISION
}

const char * ECHMET_CC errorToString(const RetCode tRet) noexcept
{
	switch (tRet) {
		ERROR_CODE_CASE(OK);
		ERROR_CODE_CASE(E_NO_MEMORY);
		ERROR_CODE_CASE(E_INVALID_ARGUMENT);
		ERROR_CODE_CASE(E_BAD_INPUT);
		ERROR_CODE_CASE(E_DATA_TOO_LARGE);
		ERROR_CODE_CASE(E_NOT_IMPLEMENTED);
		ERROR_CODE_CASE(E_LOGIC_ERROR);
		ERROR_CODE_CASE(E_NOT_FOUND);
		ERROR_CODE_CASE(E_INVALID_CONSTITUENT);
		ERROR_CODE_CASE(E_DUPLICIT_CONSTITUENTS);
		ERROR_CODE_CASE(E_INVALID_COMPLEXATION);
		ERROR_CODE_CASE(E_INVALID_COMPOSITION);
		ERROR_CODE_CASE(E_FAST_ESTIMATE_FAILURE);
		ERROR_CODE_CASE(E_NRS_FAILURE);
		ERROR_CODE_CASE(E_IS_NO_CONVERGENCE);
		ERROR_CODE_CASE(E_MISSING_PB);
		ERROR_CODE_CASE(E_NRS_NO_CONVERGENCE);
		ERROR_CODE_CASE(E_NRS_STUCK);
		ERROR_CODE_CASE(E_NRS_NO_SOLUTION);
		ERROR_CODE_CASE(E_BUFFER_CAPACITY_UNSOLVABLE);
	default:
		return "Unknown error";
	}
}

FixedStringImpl::FixedStringImpl(const std::string &stlStr) :
	m_stlStr(stlStr)
{
}

FixedStringImpl::~FixedStringImpl() noexcept
{
}

const char * ECHMET_CC FixedStringImpl::c_str() const noexcept
{
	return m_stlStr.c_str();
}

void ECHMET_CC FixedStringImpl::destroy() const noexcept
{
	delete this;
}

size_t ECHMET_CC FixedStringImpl::length() const noexcept
{
	return m_stlStr.length();
}

FixedStringImpl & ECHMET_CC FixedStringImpl::operator=(const std::string &stlStr)
{
	m_stlStr = stlStr;

	return *this;
}

bool ECHMET_CC FixedStringImpl::operator==(const FixedString &other) const noexcept
{
	try {
		const FixedStringImpl &otherImpl = dynamic_cast<const FixedStringImpl &>(other);

		return m_stlStr == otherImpl.m_stlStr;
	} catch (std::bad_cast &) {
		return false;
	}
}

bool ECHMET_CC FixedStringImpl::operator!=(const FixedString &other) const noexcept
{
	return !(*this == other);
}

FixedString::~FixedString() noexcept {}

template <> Vec<ECHMETReal>::~Vec() noexcept {}

FixedString * ECHMET_CC cpuIdentifier() ECHMET_NOEXCEPT
{
	const auto name = CPUFeatures::name();

	return createFixedString(name.c_str());
}

CPUSIMD ECHMET_CC cpuSupportedSIMD() ECHMET_NOEXCEPT
{
	const CPUFeatures::SupportedSIMD &simdInternal = CPUFeatures::SIMD();

	CPUSIMD simd;
	std::memset(&simd, 0, sizeof(CPUSIMD));

	SET_SIMD_FLAG(simd, simdInternal, SSE2);
	SET_SIMD_FLAG(simd, simdInternal, SSE3);
	SET_SIMD_FLAG(simd, simdInternal, SSSE3);
	SET_SIMD_FLAG(simd, simdInternal, SSE41);
	SET_SIMD_FLAG(simd, simdInternal, SSE42);
	SET_SIMD_FLAG(simd, simdInternal, AVX);
	SET_SIMD_FLAG(simd, simdInternal, AVX2);
	SET_SIMD_FLAG(simd, simdInternal, FMA3);

	simd.AVX512.F = simdInternal.AVX512.F;
	simd.AVX512.CD = simdInternal.AVX512.CD;
	simd.AVX512.PF = simdInternal.AVX512.PF;
	simd.AVX512.ER = simdInternal.AVX512.ER;
	simd.AVX512.VL = simdInternal.AVX512.VL;
	simd.AVX512.BW = simdInternal.AVX512.BW;
	simd.AVX512.DQ = simdInternal.AVX512.DQ;
	simd.AVX512.IFMA = simdInternal.AVX512.IFMA;
	simd.AVX512.VBM = simdInternal.AVX512.VBM;
	simd.AVX512.VBM2 = simdInternal.AVX512.VBM2;
	simd.AVX512.VNNI = simdInternal.AVX512.VNNI;
	simd.AVX512.BITALG = simdInternal.AVX512.BITALG;
	simd.AVX512.VPOPCNTDQ = simdInternal.AVX512.VPOPCNTDQ;
	simd.AVX512._4VNNIW = simdInternal.AVX512._4VNNIW;
	simd.AVX512._4FMAPS = simdInternal.AVX512._4FMAPS;
	simd.AVX512.VP2INTERSECT = simdInternal.AVX512.VP2INTERSECT;
	simd.AVX512.FP16 = simdInternal.AVX512.FP16;
	simd.AVX512.BF16 = simdInternal.AVX512.BF16;

	return simd;
}

NonidealityCorrections ECHMET_CC defaultNonidealityCorrections() noexcept
{
	return static_cast<NonidealityCorrections>(0);
}

bool ECHMET_CC nonidealityCorrectionIsSet(const NonidealityCorrections corrections, const NonidealityCorrectionsItems item) noexcept
{
	typedef typename std::underlying_type<NonidealityCorrectionsItems>::type Type;
	return corrections & static_cast<Type>(item);
}

void ECHMET_CC nonidealityCorrectionSet(NonidealityCorrections &corrections, const NonidealityCorrectionsItems item) noexcept
{
	typedef typename std::underlying_type<NonidealityCorrectionsItems>::type Type;
	corrections = corrections | static_cast<Type>(item);
}

} // namespace ECHMET
