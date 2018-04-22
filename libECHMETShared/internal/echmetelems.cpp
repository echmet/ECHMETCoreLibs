#include "../containers/echmetfixedstring_p.h"
#include "../containers/echmetvec_p.h"

#include <typeinfo>

#define _STRINGIFY(input) #input
#define ERROR_CODE_CASE(erCase) case RetCode::erCase: return _STRINGIFY(erCase)

namespace ECHMET {

RealVec * ECHMET_CC createRealVec(const size_t reserve) noexcept
{
	return createECHMETVec<ECHMETReal, false>(reserve);
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

FixedStringImpl::~FixedStringImpl()
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
