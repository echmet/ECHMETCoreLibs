#ifndef ECHMET_CAES_TOTALEQUILIBRIUM_H
#define ECHMET_CAES_TOTALEQUILIBRIUM_H

#include <echmetcaes.h>
#include <vector>
#include "funcs.h"
#ifdef ECHMET_USE_X86_EXTENSIONS
#include "vecmath/vecmath.h"
#else
#include "genericmath.h"
#endif // ECHMET_USE_X86_EXTENSIONS

namespace ECHMET {
namespace CAES {

template <typename CAESReal, InstructionSet ISet>
inline
void calculateTsAnddTsdV(CAESReal *const ECHMET_RESTRICT_PTR ts,
			 CAESReal *const ECHMET_RESTRICT_PTR dts,
			 const CAESReal &v,
			 const std::vector<CAESReal> &activityCoefficients,
			 CAESReal &X, CAESReal &dX,
			 const std::vector<CAESReal> &Ls, const size_t len,
			 const int numLow);

template <typename CAESReal, InstructionSet ISet>
inline
void calculateTsAnddTsdV(CAESReal *const ECHMET_RESTRICT_PTR ts,
			 CAESReal *const ECHMET_RESTRICT_PTR dts,
			 const CAESReal &v,
			 CAESReal &X, CAESReal &dX,
			 const std::vector<CAESReal> &Ls, const size_t len);

template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 const std::vector<double> &activityCoefficients,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len,
			 const int numLow);
template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len);

template <>
void calculateTsAnddTsdV<double, InstructionSet::FMA3>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 const std::vector<double> &activityCoefficients,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len,
			 const int numLow);
template <>
void calculateTsAnddTsdV<double, InstructionSet::FMA3>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len);

#ifndef ECHMET_DISABLE_AVX512

template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX512>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 const std::vector<double> &activityCoefficients,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len,
			 const int numLow);

template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX512>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len);

#endif // ECHMET_DISABLE_AVX512

/*!
 * Calculates total equilibrium constants from constecutive constants
 *
 * @param[in] pBs Consecutive equilibrium constants
 *
 * @return Vector of total equilibrium constants
 */
template <typename CAESReal>
std::vector<CAESReal> calculateLsBase(const std::vector<CAESReal> &pBs)
{
	std::vector<CAESReal> _Ls;
	const size_t len = pBs.size() + 1;

	_Ls.resize(len);
	_Ls[0] = 1.0;

	for (size_t idx = 1; idx < len; idx++)
		_Ls[idx] = _Ls[idx - 1] / X10(pBs[idx - 1] - 3.0); /* Use 1/L values so we can do vPow * L instead of vPow / L later on */

	return _Ls;
}

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class TotalEquilibrium;

template <typename CAESReal, InstructionSet ISet>
class TotalEquilibrium<CAESReal, ISet, true> {
public:
	typedef std::tuple<std::vector<CAESReal>, std::vector<CAESReal>> DDPack;

	TotalEquilibrium();

	/*!
	 * TotalEquilibrium c-tor
	 *
	 * @param[in] numLow Lowest equilibrium index
	 * @param[in] numHigh Highest equilibirium index
	 * @param[in] pBs Vector of consecutive equilibrium constants
	 * @param[in] concentration Analytical concentration of the constituent
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<ECHMETReal> &pBs, const size_t concentrationIndex) :
		concentrationIndex(concentrationIndex),
		Ls(calculateLs(pBs)),
		numLow(numLow),
		numHigh(numHigh)
	{
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex = SIZE_MAX);
	TotalEquilibrium(const TotalEquilibrium &other);
	TotalEquilibrium(TotalEquilibrium &&other);
	std::vector<CAESReal> concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations) const;
	DDPack distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const;
	std::vector<CAESReal> distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const;
	std::vector<CAESReal> dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const;
	std::vector<CAESReal> Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const;
	std::vector<CAESReal> Ts(const CAESReal &v, CAESReal &X) const;
	DDPack TsAnddTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X, CAESReal &dX);
	DDPack TsAnddTsdV(const CAESReal &v, CAESReal &X, CAESReal &dX);

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */
	const size_t len;

private:
	/*!
	 * Calculates total equilibrium constants from constecutive constants
	 *
	 * @param[in] pBs Consecutive equilibrium constants
	 *
	 * @return Vector of total equilibrium constants
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	std::vector<CAESReal> calculateLs(const std::vector<ECHMETReal> &pBs)
	{
		std::vector<CAESReal> _pBs;
		_pBs.reserve(pBs.size());

		for (const ECHMETReal &d : pBs)
			_pBs.emplace_back(d);

		return calculateLsBase<CAESReal>(_pBs);
	}

	std::vector<CAESReal> calculateLs(const std::vector<CAESReal> &pBs)
	{
		return calculateLsBase<CAESReal>(pBs);
	}
};

template <typename CAESReal, InstructionSet ISet>
class TotalEquilibrium<CAESReal, ISet, false> {
public:
	typedef std::tuple<const std::vector<CAESReal> &, const std::vector<CAESReal> &> DDPack;

	TotalEquilibrium();

	/*!
	 * TotalEquilibrium c-tor
	 *
	 * @param[in] numLow Lowest equilibrium index
	 * @param[in] numHigh Highest equilibirium index
	 * @param[in] pBs Vector of consecutive equilibrium constants
	 * @param[in] concentration Analytical concentration of the constituent
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<ECHMETReal> &pBs, const size_t concentrationIndex) :
		concentrationIndex(concentrationIndex),
		Ls(calculateLs(pBs)),
		numLow(numLow),
		numHigh(numHigh)
	{
		m_concentrations.resize(Ls.size());
		m_distribution.resize(Ls.size());
		m_dDistdV.resize(Ls.size());
		m_dTsdV.resize(Ls.size());
		m_Ts.resize(Ls.size());
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex = SIZE_MAX);
	TotalEquilibrium(const TotalEquilibrium &other);
	TotalEquilibrium(TotalEquilibrium &&other);
	const std::vector<CAESReal> & concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations);
	DDPack distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients);
	const std::vector<CAESReal> & distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients);
	std::vector<CAESReal> & dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X);
	const std::vector<CAESReal> & Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X);
	const std::vector<CAESReal> & Ts(const CAESReal &v, CAESReal &X);
	DDPack TsAnddTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X, CAESReal &dX);
	DDPack TsAnddTsdV(const CAESReal &v, CAESReal &X, CAESReal &dX);

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */
	const size_t len;

private:
	/*!
	 * Calculates total equilibrium constants from constecutive constants
	 *
	 * @param[in] pBs Consecutive equilibrium constants
	 *
	 * @return Vector of total equilibrium constants
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	std::vector<CAESReal> calculateLs(const std::vector<ECHMETReal> &pBs)
	{
		std::vector<CAESReal> _pBs;
		_pBs.reserve(pBs.size());

		for (const ECHMETReal &d : pBs)
			_pBs.emplace_back(d);

		return calculateLs<CAESReal>(_pBs);
	}

	std::vector<CAESReal> calculateLs(const std::vector<CAESReal> &pBs)
	{
		return calculateLsBase<CAESReal>(pBs);
	}


	std::vector<CAESReal> m_concentrations;
	std::vector<CAESReal> m_distribution;
	std::vector<CAESReal> m_dDistdV;
	std::vector<CAESReal> m_dTsdV;
	std::vector<CAESReal> m_Ts;
};

}
}

#include "totalequilibrium.hpp"

#endif // ECHMET_CAES_TOTALEQUILIBRIUM_H
