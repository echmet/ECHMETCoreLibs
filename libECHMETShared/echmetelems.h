#ifndef ECHMET_ELEMS_H
#define ECHMET_ELEMS_H

#include <echmetmodule.h>
#include <cstddef>

#ifdef ECHMET_USE_HIGH_PRECISION
#include "internal/mpreal.h"
#endif // ECHMET_USE_HIGH_PRECISION


#if __cplusplus >= 201103L
	#define _STRINGIFY_TYPENAME(theType) #theType
	#ifdef ECHMET_USE_HIGH_PRECISION
		#define IS_POD(theType) static_assert(true, "");
	#else
		#define IS_POD(theType) static_assert(std::is_pod<theType>::value, "Type " _STRINGIFY_TYPENAME(theType) " is not a POD type");
	#endif // ECHMET_USE_HIGH_PRECISION

	#ifdef ECHMET_DEBUG_OUTPUT
	#define ECHMET_DEBUG_CODE(x) { auto __df = [&](){ x; }; __df(); }
	#else
	#define ECHMET_DEBUG_CODE(x) ;
	#endif

	#include <type_traits>
	#include <cstdint>
#else
	/* Disable static checks and debugging on C++11-unaware compilers as they are not supported
	 * by the older language standards */
	#define IS_POD(theType)
	#define ECHMET_DEBUG_CODE(x)
	/* cstdint is part of the standard since C++11 */
	#include <stdint.h>

	#ifndef SIZE_MAX
	#define SIZE_MAX -1
	#endif // SIZE_MAX
#endif // C++11


/*!
 * The enclosing namespace for all functionality provided by
 * the collection of ECHMET packages
 */
namespace ECHMET {

#ifdef ECHMET_USE_HIGH_PRECISION
typedef mpfr::mpreal ECHMETReal;
#else
typedef double ECHMETReal;
#endif // ECHMET_USE_HIGH_PRECISION

/*!
 * Possible return values from public API calls.
 */
ECHMET_ST_ENUM(RetCode) {
	/* Common error codes */
	OK = 0,
	E_NO_MEMORY = 0x1,			/*!< Insufficient memory to complete operation. */
	E_INVALID_ARGUMENT = 0x2,		/*!< Argument passed to a function was invalid. */
	E_BAD_INPUT = 0x3,			/*!< Input data are malformed or nonsensical. */
	E_DATA_TOO_LARGE = 0x4,			/*!< Input data are too large to process. */
	E_NOT_IMPLEMENTED = 0x5,		/*!< Requested functionality is not implemented. */
	E_LOGIC_ERROR = 0x6,			/*!< A logic invariant was violated during exectution. */
	E_NOT_FOUND = 0x7,			/*!< Requested data was not found. */
	/* SysComp error codes */
	E_INVALID_CONSTITUENT = 0x100,		/*!< Properties of a constituent are invalid or nonsensical. */
	E_INVALID_COMPLEXATION = 0x101,		/*!< Complexation scheme is invalid or nonsensical. */
	E_DUPLICIT_CONSTITUENTS = 0x102,	/*!< There are two or more constituents of the same name. */
	/* CAES error codes */
	E_INVALID_COMPOSITION = 0x200,		/*!< System composition is ill-defined. */
	E_SOLVER_NOT_INITIALIZED = 0x201,	/*!< Solver was not initialized prior to calling <tt>solve()</tt>. */
	E_NRS_FAILURE = 0x202,			/*!< Numerical error occured in the Newton-Raphson solver. */
	E_IS_NO_CONVERGENCE = 0x203,		/*!< Solver failed to find a solution corrected for the ionic strength. */
	E_MISSING_PB = 0x204,			/*!< Complexation constant for a mixed complex form has not been entered. */
	E_NRS_NO_CONVERGENCE = 0x205,		/*!< Newton-Raphson solver failed to converge within the given number of iterations. */
	E_NRS_STUCK = 0x206,			/*!< Greatest change of X-value calculated by the Newton-Raphson solver is below the precision threshold. */
	E_NRS_NO_SOLUTION = 0x207,		/*!< System appears to have no solution. */
	E_BUFFER_CAPACITY_UNSOLVABLE = 0x208	/*!< Buffer capacity for the given system cannot be calculated. */
	ENUM_FORCE_INT32_SIZE(ECHMETRetCode)
};

/*!
 * Immutable string
 */
class FixedString {
public:
	/*!
	 * Returns CString representation of the \p FixedString.
	 *
	 * @return Pointer to the CString.
	 */
	virtual const char * ECHMET_CC c_str() const ECHMET_NOEXCEPT = 0;

	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*
	 * Returns length of the string in bytes.
	 *
	 * @return Length of the string in bytes.
	 */
	virtual size_t ECHMET_CC length() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Compares two \p FixedString s
	 *
	 * @retval true Strings have the same content.
	 * @retval false Strings have different content.
	 */
	virtual bool ECHMET_CC operator==(const FixedString &other) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Compares two \p FixedString s.
	 *
	 * @retval true Strings have different content.
	 * @retval false Strings have the same content.
	 */
	virtual bool ECHMET_CC operator!=(const FixedString &other) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~FixedString() ECHMET_NOEXCEPT = 0;
};

template <typename T>
class SKMapIterator {
public:
	/*!
	 * Advances the iterator by a given number of steps.
	 *
	 * @param[in] step Number of steps to advance the iterator by.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC advance(const int32_t step) ECHMET_NOEXCEPT = 0;

	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns the distance between two items pointed to by two iterator.
	 * The behavior is undefined if the iterators do not correspond to the same map.
	 *
	 * @return Distance between two items.
	 */
	virtual int32_t ECHMET_CC distance(const SKMapIterator<T> *other) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Fast-forwards the iterator to the last item.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC fastForward() ECHMET_NOEXCEPT = 0;

	/*!
	 * Checks if the iterator points to the last item.
	 *
	 * @retval true Iterator points to the last item.
	 * @retval false Iterator does not point to the last item.
	 */
	virtual bool ECHMET_CC hasNext() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns the key for the item the iterator points to.
	 *
	 * @return Current key.
	 */
	virtual const char * ECHMET_CC key() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Advances the iterator one step forward.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC next() ECHMET_NOEXCEPT = 0;

	/*!
	 * Rewinds the iterator one step back.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC prev() ECHMET_NOEXCEPT = 0;

	/*!
	 * Rewinds the iterator to the fisrt item.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC rewind() ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns the value for the item the iterator points to.
	 *
	 * @return Current value.
	 */
	virtual const T & ECHMET_CC value() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Advances the iterator one step forward.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC operator++() ECHMET_NOEXCEPT = 0;

	/*!
	 * Rewinds the iterator one step back.
	 *
	 * @return Pointer to self.
	 */
	virtual SKMapIterator<T> * ECHMET_CC operator--() ECHMET_NOEXCEPT = 0;

protected:
	virtual ~SKMapIterator() ECHMET_NOEXCEPT  = 0;
};

/*!
 * Template immutable Cstring-keyed map
 */
template <typename T>
class SKMap {
public:
	typedef SKMapIterator<T> Iterator;

	/*!
	 * Gets the value stored under the given key.
	 *
	 * @param[out] item Reference to variable to hold the result.
	 * @param[in] key Key in the map as CString.
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_NOT_FOUND The map does not contain the given key.
	 */
	virtual RetCode ECHMET_CC at(T &item, const char *key) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Creates a new iterator object pointing to the beginning of the map
	 *
	 * @return Iterator object for the given map
	 */
	virtual SKMapIterator<T> * ECHMET_CC begin() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Checks if the map contains the given key.
	 *
	 * @param[in] key The key to check as CString.
	 *
	 * @retval true The key is contained.
	 * @retval false The key is not contained.
	 */
	virtual bool ECHMET_CC contains(const char *key) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Creates a new iterator object pointing to the beginning of the map
	 *
	 * @return Iterator object for the given map
	 */
	virtual SKMapIterator<T> * ECHMET_CC end() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a const reference to item referenced by a given key.
	 * If the map does not contain the given key the behavior is undefined.
	 *
	 * @param[in] key Key of the requested item.
	 */
	virtual const T & ECHMET_CC operator[](const char *key) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~SKMap() ECHMET_NOEXCEPT = 0;
};

/*!
 * Template mutable Cstring-keyed map
 */
template <typename T>
class MutSKMap {
public:
	typedef SKMapIterator<T> Iterator;

	/*!
	 * Gets the value stored under the given key.
	 *
	 * @param[out] item Reference to variable to hold the result.
	 * @param[in] key Key in the map as CString.
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_NOT_FOUND The map does not contain the given key.
	 */
	virtual RetCode ECHMET_CC at(T &item, const char *key) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Creates a new iterator object pointing to the beginning of the map
	 *
	 * @return Iterator object for the given map
	 */
	virtual SKMapIterator<T> * ECHMET_CC begin() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Checks if the map contains the given key.
	 *
	 * @param[in] key The key to check as CString.
	 *
	 * @retval true The key is contained.
	 * @retval false The key is not contained.
	 */
	virtual bool ECHMET_CC contains(const char *key) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Creates a new iterator object pointing to the beginning of the map
	 *
	 * @return Iterator object for the given map
	 */
	virtual SKMapIterator<T> * ECHMET_CC end() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a reference to item referenced by a given key.
	 * If the map does not contain the given key the behavior is undefined.
	 *
	 * @param[in] key Key of the requested item.
	 */
	virtual T & ECHMET_CC item(const char *key) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a const reference to item referenced by a given key.
	 * If the map does not contain the given key the behavior is undefined.
	 *
	 * @param[in] key Key of the requested item.
	 */
	virtual const T & ECHMET_CC operator[](const char *key) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~MutSKMap() ECHMET_NOEXCEPT = 0;
};

/*!
 * Template vector
 */
template <typename T>
class Vec {
public:
	/*!
	 * Appends a vector to the vector. Vectors must be of the same type.
	 *
	 * @param[in] vec Vector to be appended.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to append the vector.
	 * @retval RetCode::E_INVALID_ARGUMENT Incompatible vector types.
	 */
	virtual RetCode ECHMET_CC append_vector(const Vec<T> *vec) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a const reference to item at a given index.
	 *
	 * @param[in] idx Item index.
	 *
	 * @return Const reference to the item.
	 */
	virtual const T & ECHMET_CC at(const size_t idx) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a const reference to the last item in the vector.
	 *
	 * @return Const reference to the item.
	 */
	virtual const T & ECHMET_CC back() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Makes a copy of the vector.
	 * If the contents of the vector are not trivially copyable
	 * the duplicated vector contains shallow copies of the items
	 * from the source vector.
	 *
	 * @return Pointer to the copy of the vector, NULL of the operation fails.
	 */
	virtual Vec<T> * ECHMET_CC duplicate() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a const reference to the first item in the vector.
	 *
	 * @return Const reference to the item.
	 */
	virtual const T & ECHMET_CC front() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Removes the last item from the vector.
	 */
	virtual void ECHMET_CC pop_back() ECHMET_NOEXCEPT = 0;

	/*!
	 * Appends an item to the vector,
	 *
	 * @param[in] item The item to be appened.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to append the item.
	 */
	virtual RetCode ECHMET_CC push_back(const T &item) ECHMET_NOEXCEPT = 0;

	/*!
	 * Prepends an item to the vector,
	 *
	 * @param[in] item The item to be prepended.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to append the item.
	 */
	virtual RetCode ECHMET_CC push_front(const T &item) ECHMET_NOEXCEPT = 0;

	/*!
	 * Releases resources claimed by the items in the vector
	 * by calling <tt>destroy()</tt> on each object.
	 * If the contents of the vector do not export such a method,
	 * calling this function is a no-op.
	 */
	virtual void ECHMET_CC releaseContents() ECHMET_NOEXCEPT = 0;

	/*!
	 * Reserves memory for the given number of items.
	 *
	 * @param[in] size Number of items to reserve memory for.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to do the reservation.
	 */
	virtual RetCode ECHMET_CC reserve(const size_t size) ECHMET_NOEXCEPT = 0;

	/*!
	 * Resizes the vector to a given size.
	 *
	 * @param[in] size The size to resize the vector to.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to resize the vector.
	 */
	virtual RetCode ECHMET_CC resize(const size_t size) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns the current number of items in the vector.
	 *
	 * @return Number of items.
	 */
	virtual size_t ECHMET_CC size() const ECHMET_NOEXCEPT = 0;

	/*!
	 * Creates a new vector which is a slice of the current vector.
	 * If the contents of the vector are not trivially copyable
	 * the duplicated vector contains shallow copies of the items
	 * from the source vector.
	 *
	 * @param[in] from First item of the slice.
	 * @param[in] to Last item of the slice. <tt>SIZE_MAX</tt> is an alias for the last item in the source vector.
	 *
	 * @return Slice of the vector, <tt>NULL</tt> if the operation failed.
	 */
	virtual Vec<T> * ECHMET_CC slice(const size_t from, size_t to = SIZE_MAX) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns a reference to the item at a given index.
	 *
	 * @param[in] idx Index of the item.
	 *
	 * @return Reference to item.
	 */
	virtual T & ECHMET_CC operator[](const size_t idx) ECHMET_NOEXCEPT = 0;

protected:
	virtual ~Vec() ECHMET_NOEXCEPT = 0;
};

typedef Vec<ECHMETReal> RealVec;

ECHMET_ST_ENUM(NonidealityCorrectionsItems) {
	CORR_DEBYE_HUCKEL = 0x1,	/* Debye-Hückel correction of stablity constants. Requires ionic strength
					   as input parameter. */
	CORR_ONSAGER_FUOSS = 0x2,	/* Onsager-Fuoss correction of ionic mobilities. Requires ionic strength
					   as input parameter. */
	CORR_VISCOSITY = 0x4		/* Viscosity correction of ionic mobilities. Requires viscosity coefficients
					   as input parameter. */
};

typedef ECHMET_ST_ENUM_UTYPE(NonidealityCorrectionsItems) NonidealityCorrections;

extern "C" {

ECHMET_API bool ECHMET_CC nonidealityCorrectionIsSet(const NonidealityCorrections corrections, const NonidealityCorrectionsItems item) ECHMET_NOEXCEPT;
ECHMET_API void ECHMET_CC nonidealityCorrectionSet(NonidealityCorrections &corrections, const NonidealityCorrectionsItems item) ECHMET_NOEXCEPT;

/*!
 * Creates an <tt>ECHMET::Vec</tt> of doubles.
 *
 * @return Pointer to the \p Vec interface, \p NULL on failure.
 */
ECHMET_API RealVec * ECHMET_CC createRealVec(const size_t reserve) ECHMET_NOEXCEPT;

/*!
 * Creates na <tt>ECHMET::FixedString</tt>.
 *
 * @return Pointer to the \p FixedString interface, \p NULL on failure.
 */
ECHMET_API FixedString * ECHMET_CC createFixedString(const char *str) ECHMET_NOEXCEPT;

/*!
 * Safe convertor from ECHMETReal internal datatype to IEEE754 \p double.
 *
 * @param[in] real ECHMETReal variable
 *
 * @return real a \p double
 */
ECHMET_API double ECHMET_CC ECHMETRealToDouble(const ECHMETReal &real) ECHMET_NOEXCEPT;

/*!
 * Returns a string representation of \p ECHMET\::RetCode return code.
 *
 * @return String representation of the return code, <tt>"Unknown error"</tt> of the value is not recognized.
 */
ECHMET_API const char * ECHMET_CC errorToString(const RetCode tRet) ECHMET_NOEXCEPT;

} // extern "C"

} // namespace ECHMET

#endif // ECHMET_ELEMS_H
