#ifndef ECHMET_VEC_P_H
#define ECHMET_VEC_P_H

#ifdef USE_ECHMET_CONTAINERS

#include <stdexcept>
#include <type_traits>
#include <vector>
#include <echmetelems.h>

namespace ECHMET {

/*!
 * Template wrapper around \p std\::vector
 */
template <typename T, bool ContentReleasable>
class VecImpl : public Vec<T> {
private:
	template<bool B>
	struct bool_ {};

public:
	typedef std::vector<T> STLVec;

	VecImpl(const std::vector<T> &stlVec = std::vector<T>()) :
		m_stlVec(stlVec)
	{
	}

	VecImpl(const std::vector<T> &&stlVec) :
		m_stlVec(stlVec)
	{
	}

	virtual ~VecImpl() noexcept override
	{
	}

	virtual RetCode ECHMET_CC append_vector(const Vec<T> *vec) noexcept override
	{
		const VecImpl<T, ContentReleasable> *vecImpl = dynamic_cast<const VecImpl<T, ContentReleasable> *>(vec);
		if (vecImpl == nullptr)
			return RetCode::E_INVALID_ARGUMENT;

		try {
			m_stlVec.insert(m_stlVec.end(), vecImpl->STL().begin(), vecImpl->STL().end());
		} catch (std::bad_alloc &) {
			return RetCode::E_NO_MEMORY;
		}

		return RetCode::OK;
	}

	virtual const T & ECHMET_CC at(const size_t idx) const noexcept override
	{
		try {
			return m_stlVec.at(idx);
		} catch (std::out_of_range &) {
			std::terminate();
		}
	}

	virtual const T & ECHMET_CC back() const noexcept override
	{
		if (m_stlVec.size() < 1)
			std::terminate();

		return m_stlVec.back();
	}

	virtual Vec<T> * ECHMET_CC duplicate() const noexcept override
	{
		return new (std::nothrow) VecImpl(m_stlVec);
	}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	virtual const T & ECHMET_CC elem(const size_t idx) const noexcept override
	{
		return m_stlVec[idx];
	}

	virtual const T & ECHMET_CC front() const noexcept override
	{
		if (m_stlVec.size() < 1)
			std::terminate();

		return m_stlVec.front();
	}

	virtual void ECHMET_CC pop_back() noexcept override
	{
		m_stlVec.pop_back();
	}

	virtual RetCode ECHMET_CC push_back(const T &item) noexcept override
	{
		try {
			m_stlVec.push_back(item);
		} catch (std::bad_alloc &) {
			return RetCode::E_NO_MEMORY;
		}

		return RetCode::OK;
	}

	virtual RetCode ECHMET_CC push_front(const T &item) noexcept override
	{
		try {
			m_stlVec.emplace(m_stlVec.begin(), item);
		} catch (std::bad_alloc &) {
			return RetCode::E_NO_MEMORY;
		}

		return RetCode::OK;
	}

	/*!
	 * If a vector contains ECHMETObject objects, calling
	 * this function will dispose of these objects. Calling this on
	 * a vector that does not contain ECHMETObject objects is a no-op.
	 */
	virtual void ECHMET_CC releaseContents() noexcept override
	{
		releaseContentsInternal(bool_<ContentReleasable>());
	}

	virtual RetCode ECHMET_CC reserve(const size_t size) noexcept override
	{
		try {
			m_stlVec.reserve(size);
		} catch (std::bad_alloc &) {
			return RetCode::E_NO_MEMORY;
		}

		return RetCode::OK;
	}

	virtual RetCode ECHMET_CC resize(const size_t size) noexcept override
	{
		try {
			m_stlVec.resize(size);
		} catch (std::bad_alloc &) {
			return RetCode::E_NO_MEMORY;
		}

		return RetCode::OK;
	}

	virtual size_t ECHMET_CC size() const noexcept override
	{
		return m_stlVec.size();
	}

	virtual Vec<T> * ECHMET_CC slice(const size_t from, size_t to = SIZE_MAX) const noexcept override
	{
		if (to == SIZE_MAX)
			to = m_stlVec.size() - 1;

		if (from > to)
			return nullptr;
		if (to >= m_stlVec.size())
			return nullptr;


		std::vector<T> target;

		try {
			target.resize(to - from + 1);
			std::copy(m_stlVec.cbegin() + from, m_stlVec.cbegin() + to, target.begin());

			return new VecImpl(std::move(target));
		} catch (std::bad_alloc &) {
			return nullptr;
		}
	}

	virtual T & ECHMET_CC operator[](const size_t idx) noexcept override
	{
		return m_stlVec[idx];
	}

	std::vector<T> & STL() noexcept
	{
		return m_stlVec;
	}

	const std::vector<T> & STL() const noexcept
	{
		return m_stlVec;
	}

private:
	void contentReleaser(T object, bool_<true>) noexcept
	{
		object->destroy();
	}

	void releaseContentsInternal(bool_<true>) noexcept
	{
		for (size_t idx = 0; idx < m_stlVec.size(); idx++) {
			T object = m_stlVec[idx];

			contentReleaser(object, bool_<ContentReleasable>());
		}
	}

	void releaseContentsInternal(bool_<false>) noexcept
	{
		return;
	}

	std::vector<T> m_stlVec;
};

template <typename T, bool ContentReleasable>
VecImpl<T, ContentReleasable> * ECHMET_CC createECHMETVec(const size_t reserve) noexcept
{
	try {
		VecImpl<T, ContentReleasable> *v = new VecImpl<T, ContentReleasable>();

		if (reserve > 0)
			v->reserve(reserve);

		return v;
	} catch (std::bad_alloc &) {
		return nullptr;
	}

}

template <typename T>
Vec<T>::~Vec() noexcept {}

} // namespace ECHMET

#endif// USE_ECHMET_CONTAINERS

#endif // ECHMET_VEC_P_H
