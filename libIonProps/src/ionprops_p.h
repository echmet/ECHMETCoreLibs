#ifndef ECHMET_IONPROPS_IONPROPS_P_H
#define ECHMET_IONPROPS_IONPROPS_P_H

#include <echmetionprops.h>
#include <cstddef>
#include <vector>

#include <echmetmodule.h>
#include <internal/mpreal.h>

namespace ECHMET {
namespace IonProps {

template <bool>
struct IsNative {};

template <typename IPReal>
static double IPRealToDouble(const IPReal &v);

template <>
double IPRealToDouble(const mpfr::mpreal &v)
{
	return v.toDouble();
}

template <>
double IPRealToDouble(const double &v)
{
	return v;
}

template <bool>
struct IsSame {};

template <typename From, typename To>
static To __IPRealToECHMETReal(const From &f, IsSame<true>)
{
	return f;
}

template <typename From, typename To>
static To __IPRealToECHMETReal(const From &f, IsSame<false>);

template <>
double __IPRealToECHMETReal<mpfr::mpreal, double>(const mpfr::mpreal &f, IsSame<false>)
{
	return f.toDouble();
}

template <>
mpfr::mpreal __IPRealToECHMETReal<double, mpfr::mpreal>(const double &f, IsSame<false>)
{
	return mpfr::mpreal(f);
}

template <typename From, typename To = ECHMETReal>
static To IPRealToECHMETReal(const From &f)
{
	return __IPRealToECHMETReal<From, To>(f, IsSame<std::is_same<From, To>::value>());
}

/*!
 * Dynamically allocated 2D array flattened into 1D array.
 */
template <typename T>
class Array2D {
public:
	Array2D(const size_t rows, const size_t columns) :
		m_length(rows * columns),
		m_columns(columns)
	{
		m_array = new T[m_length];
	}

	~Array2D()
	{
		delete[] m_array;
	}

	const T & at(const size_t row, const size_t column) const
	{
		return m_array[index(row, column)];
	}

	T & operator()(const size_t row, const size_t column)
	{
		return m_array[index(row, column)];
	}

private:
	T *m_array;
	const size_t m_length;
	const size_t m_columns;

	inline size_t index(const size_t row, const size_t column) const
	{
		return m_columns * row + column;
	}
};

/*!
 * Internal representation of a charged chemical species.
 */
template <typename IPReal>
class Ion {
public:
	/*!
	 * \p Ion c-tor.
	 *
	 * @param[in] ionicMobilityIndex Index in the \p ChemicalSystem\::ionicMobilities vector where the ionic mobility of the ion is stored.
	 * @param[in] concentration Ionic concentration of the ion in the system in <tt>mmol/dm<sup>3</sup></tt>.
	 * @param[in] charge Electric charge of the ion.
	 * @param[in] limitMobility Limit electrophoretic mobility of the ion in <tt>V.s/m<sup>2</sup> \. 10e9</tt>.
	 */
	Ion(const size_t ionicMobilityIndex, const IPReal &concentration, const int charge, const IPReal &limitMobility) :
		ionicMobilityIndex{ionicMobilityIndex},
		concentration{concentration},
		charge{charge},
		limitMobility{limitMobility * 1.0e-9}
	{}

	const size_t ionicMobilityIndex;	/*<! Index in the \p ChemicalSystem\::ionicMobilities vector where the ionic mobility of the ion is stored. */
	const IPReal concentration;		/*<! Ionic concentration of the ion in the system in <tt>mmol/dm<sup>3</sup></tt>. */
	const int charge;			/*<! Electric charge of the ion. */
	const IPReal limitMobility;		/*<! Limit electrophoretic mobulity of the ion in <tt>V.s/m<sup>2</sup></tt>. */
};

/*!
 * Internal implementation of the comutation context
 */
template <typename IPReal>
class ComputationContextImpl : public ComputationContext {
public:
	/* \p ComputationContextImpl c-tor.
	 *
	 * @param[in] chemSystem The corresponding chemical system.
	 * @param[in] analyticalConcentrations Poitner to the vector of analytical concentrations of all compounds in the system.
	 * @param[in,out] calcProps Corresponding \p SysComp\::CalculatedProperties solved by \p CAES.
	 */
	explicit ComputationContextImpl(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations,
					SysComp::CalculatedProperties &calcProps) :
		chemSystem{chemSystem},
		analyticalConcentrations{analyticalConcentrations},
		calcProps{calcProps},
		ionicConcentrations{std::vector<IPReal>{}}
	{}

	/* \p ComputationContextImpl c-tor.
	 * This overload is available only when the \p ECHMETReal is typedefed to an IEEE754 floating point type.
	 *
	 * @param[in] Vector of ionic concentrations represented as \p IPReal type instead of \p ECHMETReal type.
	 * @param[in] chemSystem The corresponding chemical system.
	 * @param[in] analyticalConcentrations Poitner to the vector of analytical concentrations.
	 * @param[in,out] calcProps \p SysComp\::CalculatedProperties solved by \p CAES.
	 */
	template <typename E = ECHMETReal, typename = typename std::enable_if<std::is_floating_point<E>::value>::type>
	explicit ComputationContextImpl(const std::vector<IPReal> &ionicConcentrations, const SysComp::ChemicalSystem &chemSystem,
					const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) :
		chemSystem{chemSystem},
		analyticalConcentrations{analyticalConcentrations},
		calcProps{calcProps},
		ionicConcentrations{ionicConcentrations}
	{}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	const SysComp::ChemicalSystem &chemSystem;		/*!< The corresponding chemical system. */
	const RealVec *analyticalConcentrations;		/*!< Pointer to vector of analytical concentrations of all compounds in the system. */
	SysComp::CalculatedProperties &calcProps;		/*!< Corresponding \p SysComp\::CalculatedProperties solved by \p CAES. */
	const std::vector<IPReal> ionicConcentrations;		/*!< Vector of ionic concentrations represented as \p IPReal. This vector is empty unless
								     the object was created using the overloaded c-tor. */
};

template <typename T>
inline constexpr T cxpow(const T &x, const unsigned int exp)
{
	return (exp > 0) ? (x * (cxpow(x, exp - 1))) : 1;
}

template <typename T>
int sgn(const T &val)
{
	return (T(0) < val) - (val < T(0));
}

} // namespace IonProps
} // namespace ECHMET


#include "ionprops.hpp"

#endif // ECHMET_IONPROPS_IONPROPS_P_H
