#ifndef ECHMET_FIXEDSTRING_P_H
#define ECHMET_FIXEDSTRING_P_H

#ifdef USE_ECHMET_CONTAINERS

#include <string>
#include <echmetelems.h>

namespace ECHMET {

/*!
 * Wrapper around \p std\::string
 */
class FixedStringImpl : public FixedString {
public:
	FixedStringImpl(const std::string &stlStr = std::string());
	virtual ~FixedStringImpl() noexcept override;
	virtual const char * ECHMET_CC c_str() const noexcept override;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual size_t ECHMET_CC length() const noexcept override;

	virtual bool ECHMET_CC operator==(const FixedString &other) const noexcept override;
	virtual bool ECHMET_CC operator!=(const FixedString &other) const noexcept override;

	FixedStringImpl & ECHMET_CC operator=(const std::string &stlStr);
	bool ECHMET_CC operator==(const std::string &stlstr);

private:
	std::string m_stlStr;
};

} //namespace ECHMET

#endif // USE_ECHMET_CONTAINERS

#endif // ECHMET_FIXEDSTRING_P_H
