#ifndef ECHMET_MATH_H
#define ECHMET_MATH_H

#include <echmetelems.h>
#include <cmath>

#include "mpreal.h"

namespace ECHMET {

/*!
 * Internal template math functions
 */
namespace VMath {

template <typename XReal>
static XReal abs(const XReal &arg);

template <typename XReal>
static bool isinf(const XReal &arg);

template <typename XReal>
static bool isnan(const XReal &arg);

template <typename XReal>
static XReal log10(const XReal &arg);

template <typename XReal, typename Base, typename Exponent>
static XReal pow(const Base &base, const Exponent &exponent);

template <typename XReal>
static XReal sqrt(const XReal &arg);

/* ABS */
template <>
mpfr::mpreal abs(const mpfr::mpreal &arg)
{
	return mpfr::abs(arg);
}

template <>
double abs(const double &arg)
{
	return ::std::abs(arg);
}

template <>
int abs(const int &arg)
{
	return ::std::abs(arg);
}
/**/

/* ISINF */
template <>
bool isinf(const mpfr::mpreal &arg)
{
	return mpfr::isinf(arg);
}

template <>
bool isinf(const double &arg)
{
	return ::std::isinf(arg);
}
/**/

/* ISNAN */
template <>
bool isnan(const mpfr::mpreal &arg)
{
	return mpfr::isnan(arg);
}

template <>
bool isnan(const double &arg)
{
	return ::std::isnan(arg);
}
/**/

/* LOG10 */
template <>
mpfr::mpreal log10(const mpfr::mpreal &arg)
{
	return mpfr::log10(arg);
}

template <>
double log10(const double &arg)
{
	return ::std::log10(arg);
}
/**/

/* POW */
template <>
mpfr::mpreal pow<mpfr::mpreal, mpfr::mpreal, mpfr::mpreal>(const mpfr::mpreal &base, const mpfr::mpreal &exponent)
{
	return mpfr::pow(base, exponent);
}

template <>
mpfr::mpreal pow<mpfr::mpreal, int, mpfr::mpreal>(const int &base, const mpfr::mpreal &exponent)
{
	return mpfr::pow(base, exponent);
}

template <>
mpfr::mpreal pow<mpfr::mpreal, mpfr::mpreal, int>(const mpfr::mpreal &base, const int &exponent)
{
	return mpfr::pow(base, exponent);
}

template <>
mpfr::mpreal pow<mpfr::mpreal, double, int>(const double &base, const int &exponent)
{
	return mpfr::pow(mpfr::mpreal(base), exponent);
}

template <>
mpfr::mpreal pow<mpfr::mpreal, double, mpfr::mpreal>(const double &base, const mpfr::mpreal &exponent)
{
	return mpfr::pow(mpfr::mpreal(base), exponent);
}

template <>
mpfr::mpreal pow<mpfr::mpreal, int, int>(const int &base, const int &exponent)
{
	return mpfr::pow(base, exponent);
}

template <>
double pow<double, double, double>(const double &base, const double &exponent)
{
	return ::std::pow(base, exponent);
}

template <>
double pow<double, int, double>(const int &base, const double &exponent)
{
	return ::std::pow(base, exponent);
}

template <>
double pow<double, double, int>(const double &base, const int &exponent)
{
	return ::std::pow(base, exponent);
}

template <>
double pow<double, int, int>(const int &base, const int &exponent)
{
	return ::std::pow(base, exponent);
}

template <>
int pow<int, int, int>(const int &base, const int &exponent)
{
	return ::std::pow(base, exponent);
}

/**/

/* SQRT */
template <>
mpfr::mpreal sqrt(const mpfr::mpreal &arg)
{
	return mpfr::sqrt(arg);
}

template <>
double sqrt(const double &arg)
{
	return ::std::sqrt(arg);
}
/**/

} // namespace VMath
} // namespace ECHMET

#endif // ECHMET_MATH_H
