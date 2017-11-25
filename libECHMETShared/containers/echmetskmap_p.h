#ifndef ECHMET_SKMAP_H
#define ECHMET_SKMAP_H

#ifdef USE_ECHMET_CONTAINERS

#include <echmetelems.h>
#include <map>
#include "echmetfixedstring_p.h"

namespace ECHMET {

template <typename T>
class SKMapIteratorImpl : public SKMapIterator<T> {
public:
	SKMapIteratorImpl<T>(const std::map<std::string, T> &skMap)
	{
		m_begin = skMap.cbegin();
		m_end = skMap.cend();
		m_iterator = m_begin;
	}

	virtual ~SKMapIteratorImpl<T>() noexcept override
	{
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC advance(const int32_t step) noexcept override
	{
		std::advance(m_iterator, step);

		return this;
	}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	virtual int32_t ECHMET_CC distance(const SKMapIterator<T> *other) const noexcept
	{
		const SKMapIteratorImpl<T> *otherImpl = dynamic_cast<const SKMapIteratorImpl<T> *>(other);

		if (otherImpl == nullptr)
			return -1;

		return std::distance(m_iterator, otherImpl->m_iterator);
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC fastForward() noexcept override
	{
		m_iterator = m_end;

		return this;
	}

	virtual bool ECHMET_CC hasNext() const noexcept override
	{
		return m_iterator != m_end;
	}

	virtual const char * ECHMET_CC key() const noexcept override
	{
		if (m_iterator == m_end)
			return nullptr;

		return m_iterator->first.c_str();
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC next() noexcept override
	{
		m_iterator = std::next(m_iterator);

		return this;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC prev() noexcept override
	{
		m_iterator = std::prev(m_iterator);

		return this;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC rewind() noexcept override
	{
		m_iterator = m_begin;

		return this;
	}

	virtual const T & ECHMET_CC value() const noexcept override
	{
		if (m_iterator == m_end)
			std::terminate();

		return m_iterator->second;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC operator++() noexcept override
	{
		m_iterator++;

		return this;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC operator--() noexcept override
	{
		m_iterator--;

		return this;
	}

private:
	typename std::map<std::string, T>::const_iterator m_iterator;
	typename std::map<std::string, T>::const_iterator m_begin;
	typename std::map<std::string, T>::const_iterator m_end;
};

template <typename T>
SKMapIterator<T>::~SKMapIterator() noexcept {}

/*!
 * string-keyed template ordered map
 */
template <typename T>
class SKMapImpl : public SKMap<T> {
public:
	typedef std::map<std::string, T> STLMap;

	SKMapImpl()
	{
	}

	virtual ~SKMapImpl() noexcept override
	{
	}

	virtual RetCode ECHMET_CC at(T &item, const char *key) const noexcept override
	{
		const std::string stlKey(key);

		try {
			item = m_map.at(stlKey);

			return RetCode::OK;
		} catch (std::out_of_range &) {
			return RetCode::E_NOT_FOUND;
		}
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC begin() const noexcept override
	{
		return new (std::nothrow) SKMapIteratorImpl<T>(m_map);
	}

	virtual bool ECHMET_CC contains(const char *key) const noexcept override
	{
		const std::string &stlKey(key);

		typename std::map<std::string, T>::const_iterator it = m_map.find(stlKey);

		return it != m_map.cend();
	}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC end() const noexcept override
	{
		SKMapIteratorImpl<T> *it = new (std::nothrow) SKMapIteratorImpl<T>(m_map);
		if (it == nullptr)
			return nullptr;

		it->rewind();

		return it;
	}

	/*!
	 * Unsafe function that returns a direct reference to the object
	 * in the map
	 */
	virtual const T & ECHMET_CC operator[](const char *key) const noexcept override
	{
		return m_map.at(key);
	}

	STLMap & STL()
	{
		return m_map;
	}

	const STLMap & cSTL() const

	{
		return m_map;
	}

private:
	STLMap m_map;
};

template <typename T>
SKMap<T>::~SKMap() noexcept {}

template <typename T>
class MutSKMapImpl : public MutSKMap<T> {
public:
	typedef std::map<std::string, T> STLMap;

	MutSKMapImpl()
	{
	}

	virtual ~MutSKMapImpl() noexcept override
	{
	}

	virtual RetCode ECHMET_CC at(T &item, const char *key) const noexcept override
	{
		const std::string stlKey(key);

		try {
			item = m_map.at(stlKey);

			return RetCode::OK;
		} catch (std::out_of_range &) {
			return RetCode::E_NOT_FOUND;
		}
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC begin() const noexcept override
	{
		return new (std::nothrow) SKMapIteratorImpl<T>(m_map);
	}

	virtual bool ECHMET_CC contains(const char *key) const noexcept override
	{
		const std::string &stlKey(key);

		typename std::map<std::string, T>::const_iterator it = m_map.find(stlKey);

		return it != m_map.cend();
	}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	virtual SKMapIteratorImpl<T> * ECHMET_CC end() const noexcept override
	{
		SKMapIteratorImpl<T> *it = new (std::nothrow) SKMapIteratorImpl<T>(m_map);
		if (it == nullptr)
			return nullptr;

		it->rewind();

		return it;
	}

	/*!
	 * Unsafe function that returns a direct reference to the object
	 * in the map
	 */
	virtual T & ECHMET_CC item(const char *key) noexcept override
	{
		return m_map.at(key);
	}

	/*!
	 * Unsafe function that returns a direct reference to the object
	 * in the map
	 */
	virtual const T & ECHMET_CC operator[](const char *key) const noexcept override
	{
		return m_map.at(key);
	}

	STLMap & STL()
	{
		return m_map;
	}

	const STLMap & cSTL() const

	{
		return m_map;
	}

private:
	STLMap m_map;
};

template <typename T>
MutSKMap<T>::~MutSKMap() noexcept {}

/*!
 * Convenience function that returns
 * a pointer to SKMapImpl
 */
template <typename T>
SKMap<T> * createSKMap() noexcept
{
	return new (std::nothrow) SKMapImpl<T>();
}

/*!
 * Convenience function that returns
 * a pointer to MutSKMapImpl
 */
template <typename T>
SKMap<T> * createMutSKMap() noexcept
{
	return new (std::nothrow) MutSKMapImpl<T>();
}

} // namespace ECHMET

#endif // USE_ECHMET_CONTAINERS

#endif // ECHMET_SKMAP_H
