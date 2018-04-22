#ifndef ECHMET_CAES_TOTALEQUILIBRIUM_H
#define ECHMET_CAES_TOTALEQUILIBRIUM_H

#include <echmetcaes.h>
#include <vector>
#include "funcs.h"

namespace ECHMET {
namespace CAES {

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
	std::vector<CAESReal> _pBs;
	std::vector<CAESReal> _Ls;

	_pBs.resize(pBs.size());
	std::copy(pBs.cbegin(), pBs.cend(), _pBs.begin());

	_pBs.emplace_back(0.0);

	const size_t len = _pBs.size();

	_Ls.resize(len);

	for (size_t idx = 0; idx < len; idx++) {
		auto calcL = [](const std::vector<CAESReal> &pXs, const size_t to) {
			CAESReal L = 1.0;

			/* TODO: Optimize this!!! */
			size_t idx = pXs.size() - 1;
			for (;;) {
				L *= X10(pXs.at(idx) - 3.0);

				if (idx == to || idx == 0)
					break;
				else
					idx--;
			}

			ECHMET_DEBUG_CODE(fprintf(stderr, "L = %g\n", CAESRealToDouble(L)));

			return L;
		};

		_Ls[idx] = calcL(_pBs, idx);
	}

	return _Ls;
}

class TotalEquilibriumBase {};

/*!
 * Description of total distribution equilibria for a given constituent
 */
template <typename CAESReal, bool ThreadSafe>
class TotalEquilibrium : public TotalEquilibriumBase {
private:
};

template <typename CAESReal>
class TotalEquilibrium<CAESReal, true> : public TotalEquilibriumBase {
public:
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
	std::vector<CAESReal> concentrations(const CAESReal &v, const RealVec *analyticalConcentrations) const;
	std::vector<CAESReal> distribution(const CAESReal &v) const;
	std::vector<CAESReal> Ts(const CAESReal &v, CAESReal &X) const;

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */

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

template <typename CAESReal>
class TotalEquilibrium<CAESReal, false> : public TotalEquilibriumBase {
public:
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
		m_Ts.resize(Ls.size());
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex = SIZE_MAX);
	TotalEquilibrium(const TotalEquilibrium &other);
	TotalEquilibrium(TotalEquilibrium &&other);
	const std::vector<CAESReal> & concentrations(const CAESReal &v, const RealVec *analyticalConcentrations);
	const std::vector<CAESReal> & distribution(const CAESReal &v);
	const std::vector<CAESReal> & Ts(const CAESReal &v, CAESReal &X);

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */

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
	std::vector<CAESReal> m_Ts;
};

}
}

#include "totalequilibrium.hpp"

#endif // ECHMET_CAES_TOTALEQUILIBRIUM_H
