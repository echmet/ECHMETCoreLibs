#ifndef ECHMET_IONPROPS_IONPROPS_P_H
#define ECHMET_IONPROPS_IONPROPS_P_H

#include <echmetionprops.h>
#include <cstddef>
#include <vector>
#include <Eigen/Dense>

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
	Ion(const size_t ionicMobilityIndex, const size_t ionicConcentrationIndex, const int charge, const IPReal &limitMobility) :
		ionicMobilityIndex{ionicMobilityIndex},
		ionicConcentrationIndex{ionicConcentrationIndex},
		charge{charge},
		limitMobility{limitMobility * 1.0e-9}
	{}

	const size_t ionicMobilityIndex;	/*<! Index into the \p ChemicalSystem\::ionicMobilities vector where the ionic mobility of the ion is stored. */
	const size_t ionicConcentrationIndex;	/*<! Index into the vector of ionic concentrations vector that will be passed to the calculator function. */
	const int charge;			/*<! Electric charge of the ion. */
	const IPReal limitMobility;		/*<! Limit electrophoretic mobulity of the ion in <tt>V.s/m<sup>2</sup></tt>. */
};

/*!
 * Builds a vector of Ions. Uncharged ionic forms are not included in the vector.
 * This a templated function working with \p IPReal type.
 *
 * @param[in] ifVec Vector of \p SysComp::IonicForm objects for the corresponding chemical system.
 * @param[in] icVec Vector of ionic concentrations of all ionic forms in the system as \p IPReal.
 * @param[in] calcProps \p SysComp::CalculatedProperties for the corresponding chemical system.
 *
 * @return \p std::vector of \p Ions.
 */
template <typename IPReal>
std::vector<Ion<IPReal>> makeIonVector(const SysComp::IonicFormVec *ifVec)
{
	std::vector<Ion<IPReal>> ions;

	ions.reserve(ifVec->size());

	for (size_t idx = 0; idx < ifVec->size(); idx++) {
		const SysComp::IonicForm *iF = ifVec->at(idx);
		const ECHMETReal limitMobility = iF->limitMobility;

		if (iF->totalCharge == 0)
			continue;
#ifdef IONPROPS_DISABLE_COMPLEX_ONSFUO
		if (iF->ligand != nullptr)
			continue;
#endif // IONPROPS_DISABLE_COMPLEX_ONSFUO

		ions.emplace_back(iF->ionicMobilityIndex, iF->ionicConcentrationIndex, iF->totalCharge, limitMobility);
	}

	return ions;
}

/*!
 * Internal implementation of the comutation context
 */
template <typename IPReal>
class ComputationContextImpl : public ComputationContext {
public:
	typedef Eigen::Matrix<IPReal, Eigen::Dynamic, Eigen::Dynamic, static_cast<int>(Eigen::ColMajor) | static_cast<int>(Eigen::AutoAlign)> Matrix;

	class OnsagerFuossPack {
	public:
		OnsagerFuossPack() :
			m_releaseOnFinalize{false}
		{}

		OnsagerFuossPack(const size_t N, const bool releaseOnFinalize) :
			m_releaseOnFinalize{releaseOnFinalize}
		{
			R = Matrix(N, 6);
			H = Matrix(N, N);
			mu_I = std::vector<IPReal>(N);
		}

		void finalize()
		{
			if (m_releaseOnFinalize)
				delete this;
		}

		Matrix R;
		Matrix H;
		std::vector<IPReal> mu_I;

		OnsagerFuossPack & operator=(OnsagerFuossPack &&other) noexcept
		{
			const_cast<Matrix&>(R) = std::move(other.R);
			const_cast<Matrix&>(H) = std::move(other.H);
			const_cast<std::vector<IPReal>&>(mu_I) = std::move(other.mu_I);
			const_cast<bool&>(m_releaseOnFinalize) = other.m_releaseOnFinalize;

			return *this;
		}

	private:
		const bool m_releaseOnFinalize;
	};

	/* \p ComputationContextImpl c-tor.
	 *
	 * @param[in] chemSystem The corresponding chemical system.
	 * @param[in] analyticalConcentrations Poitner to the vector of analytical concentrations of all compounds in the system.
	 * @param[in,out] calcProps Corresponding \p SysComp\::CalculatedProperties solved by \p CAES.
	 */
	explicit ComputationContextImpl(const SysComp::ChemicalSystem &chemSystem, const Options options) :
		chemSystem{chemSystem},
		ions{makeIonVector<IPReal>(chemSystem.ionicForms)},
		ionicConcentrations{std::vector<IPReal>{}},
		m_isThreadUnsafe{(options & Options::DISABLE_THREAD_SAFETY) != 0}
	{
		prepareOnsagerFuossPack();
	}

	/* \p ComputationContextImpl c-tor.
	 * This overload is available only when the \p ECHMETReal is typedefed to an IEEE754 floating point type.
	 *
	 * @param[in] Vector of ionic concentrations represented as \p IPReal type instead of \p ECHMETReal type.
	 * @param[in] chemSystem The corresponding chemical system.
	 * @param[in] analyticalConcentrations Poitner to the vector of analytical concentrations.
	 * @param[in,out] calcProps \p SysComp\::CalculatedProperties solved by \p CAES.
	 */
	template <typename E = ECHMETReal, typename = typename std::enable_if<std::is_floating_point<E>::value>::type>
	explicit ComputationContextImpl(const std::vector<IPReal> &ionicConcentrations, const SysComp::ChemicalSystem &chemSystem, const Options options) :
		chemSystem{chemSystem},
		ions{makeIonVector<IPReal>(chemSystem.ionicForms)},
		ionicConcentrations{ionicConcentrations},
		m_isThreadUnsafe{(options & Options::DISABLE_THREAD_SAFETY) != 0}
	{
		prepareOnsagerFuossPack();
	}

	virtual void ECHMET_CC destroy() const noexcept override
	{
		delete this;
	}

	OnsagerFuossPack * onsagerFuossPack()
	{
		if (m_isThreadUnsafe)
			return &m_onsFuoPack;
		else {
			const size_t N = ions.size();

			return new OnsagerFuossPack(N, true);
		}
	}

	const SysComp::ChemicalSystem &chemSystem;		/*!< The corresponding chemical system. */
	const std::vector<Ion<IPReal>> ions;			/*!< Vector of all charged ions used by Onsager-Fuoss correction workers */
	const std::vector<IPReal> ionicConcentrations;		/*!< Vector of ionic concentrations represented as \p IPReal. This vector is empty unless
								     the object was created using the overloaded c-tor. */

private:
	void prepareOnsagerFuossPack()
	{
		if (!m_isThreadUnsafe)
			return;

		const size_t N = ions.size();

		m_onsFuoPack = OnsagerFuossPack(N, false);
	}

	const bool m_isThreadUnsafe;

	OnsagerFuossPack m_onsFuoPack;
};

template <typename T>
inline constexpr T cxpow(const T &x, const unsigned int exp)
{
	return (exp > 0) ? (x * (cxpow(x, exp - 1))) : 1;
}

template <typename T>
static
int sgn(const T &val)
{
	return (T(0) < val) - (val < T(0));
}

} // namespace IonProps
} // namespace ECHMET


#include "ionprops.hpp"

#endif // ECHMET_IONPROPS_IONPROPS_P_H
