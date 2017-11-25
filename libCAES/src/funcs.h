#ifndef ECHMET_CAES_FUNCS_H
#define ECHMET_CAES_FUNCS_H

#include "types.h"
#include <internal/echmetmath_internal.h>

#define CVI(v, idx) v(idx)

namespace ECHMET {
namespace CAES {

namespace {

template <typename CAESReal>
static double CAESRealToDouble(const CAESReal &v);

template <>
double CAESRealToDouble(const mpfr::mpreal &v)
{
	return v.toDouble();
}

template <>
double CAESRealToDouble(const double &v)
{
	return v;
}

template <bool>
struct IsSame {};

template <typename From, typename To>
static To __CAESRealToECHMETReal(const From &f, IsSame<true>)
{
	return f;
}

template <typename From, typename To>
static To __CAESRealToECHMETReal(const From &f, IsSame<false>);

template <>
double __CAESRealToECHMETReal<mpfr::mpreal, double>(const mpfr::mpreal &f, IsSame<false>)
{
	return f.toDouble();
}

template <>
mpfr::mpreal __CAESRealToECHMETReal<double, mpfr::mpreal>(const double &f, IsSame<false>)
{
	return mpfr::mpreal(f);
}

template <typename From, typename To = ECHMETReal>
static To CAESRealToECHMETReal(const From &f)
{
	return __CAESRealToECHMETReal<From, To>(f, IsSame<std::is_same<From, To>::value>());
}

}

template<typename T>
static constexpr T cxpow(const T &x, const unsigned int exp)
{
	return (exp > 0) ? (x * (cxpow(x, exp - 1))) : 1;
}

template <typename CAESReal>
static constexpr CAESReal electroneturalityPrecision();

template <>
double constexpr electroneturalityPrecision<double>()
{
	return 1.0e-13;
}

template <>
mpfr::mpreal electroneturalityPrecision<mpfr::mpreal>()
{
	return mpfr::mpreal("1.0e-60");
}

template<typename T>
static inline T pX(const T &t)
{
	return -VMath::log10(t);
}

template<typename T>
static inline constexpr int sgn(const T &val)
{
	return (T(0) < val) - (val < T(0));
}

template<typename T>
static inline T X10(const T &t)
{
	return VMath::pow<T>(10.0, -t);
}

/*!
 * Checks if a given vector of ligand ionic forms contains the given ionic form
 *
 * @param[in] clVec Vector of ligand ionic forms to inspect
 * @param[in] name Name of the ligand ionic form that is being searched for
 * @param[in, out] cl Copy of the ligand ionic form contained in the vector. This variable
 *                 is set only if the given ionic form is found in the vector.
 *
 * @return True if the ionic form is found in the vector, false otherwise
 */
template <typename CAESReal>
static bool isLigandIFContained(const ContainedLigandIonicFormVec<CAESReal> &clVec, const std::string &name, ContainedLigandIonicForm<CAESReal> &cl)
{
	for (const ContainedLigandIonicForm<CAESReal> &c : clVec) {
		if (c.lIF->name == name) {
			cl = c;
			return true;
		}
	}

	return false;
}

template<typename T>
static void releasePointerContainer(T *&container)
{
	if (container != nullptr) {
		for (const auto *i : *container)
			delete i;

		delete container;
		container = nullptr;
	}
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_FUNCS_H
