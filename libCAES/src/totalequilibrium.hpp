#ifndef ECHMET_CAES_TOTALEQUILIBRIUM_HPP
#define ECHMET_CAES_TOTALEQUILIBRIUM_HPP

namespace ECHMET {
namespace CAES {

/*!
 * TotalEquilibrium c-tor.
 *
 * @param[in] numLow Lowest equilibrium index.
 * @param[in] numHigh Highest equilibirium index.
 * @param[in] pBs Vector of consecutive equilibrium constants.
 * @param[in] concentration Analytical concentration of the constituent.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex) :
	concentrationIndex(concentrationIndex),
	Ls(calculateLs(pBs)),
	numLow(numLow),
	numHigh(numHigh)
{
}

/*!
 * TotalEquilibrium c-tor.
 *
 * @param[in] numLow Lowest equilibrium index.
 * @param[in] numHigh Highest equilibirium index.
 * @param[in] pBs Vector of consecutive equilibrium constants.
 * @param[in] concentration Analytical concentration of the constituent.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex) :
	concentrationIndex(concentrationIndex),
	Ls(calculateLs(pBs)),
	numLow(numLow),
	numHigh(numHigh)
{
	m_concentrations.resize(Ls.size());
	m_dDistdV.resize(Ls.size());
	m_distribution.resize(Ls.size());
	m_dTsdV.resize(Ls.size());
	m_Ts.resize(Ls.size());
}

template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(const TotalEquilibrium &other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh)
{
}

template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(const TotalEquilibrium &other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh),
	m_concentrations(other.m_concentrations),
	m_dDistdV(other.m_dDistdV),
	m_distribution(other.m_distribution),
	m_dTsdV(other.m_dTsdV),
	m_Ts(other.m_Ts)
{
}


/*!
 * TotalEquilibrium move c-tor.
 *
 * @param[in] other TotalEquilibrium object being moved.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(TotalEquilibrium &&other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh)
{
}

/*!
 * TotalEquilibrium move c-tor.
 *
 * @param[in] other TotalEquilibrium object being moved.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(TotalEquilibrium &&other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh),
	m_concentrations(std::move(other.m_concentrations)),
	m_dDistdV(std::move(other.m_dDistdV)),
	m_distribution(std::move(other.m_distribution)),
	m_dTsdV(std::move(other.m_dTsdV)),
	m_Ts(std::move(other.m_Ts))
{
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::concentrations(const CAESReal &v, const RealVec *analyticalConcentrations) const
{
	CAESReal X = 0.0;
	std::vector<CAESReal> concentrations{};
	const auto ts = Ts(v, X); /* Recalculate Ts */

	concentrations.resize(ts.size());
	size_t ctr = 0;
	const CAESReal &c = analyticalConcentrations->elem(concentrationIndex);
	for (const CAESReal &T : ts) {
		const CAESReal fC = c * T / X;

		concentrations[ctr] = fC;

		ctr++;
	}

	return concentrations;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::concentrations(const CAESReal &v, const RealVec *analyticalConcentrations)
{
	CAESReal X = 0.0;
	Ts(v, X); /* Recalculate Ts */

	size_t ctr = 0;
	const CAESReal &c = analyticalConcentrations->elem(concentrationIndex);
	for (const CAESReal &T : m_Ts) {
		const CAESReal fC = c * T / X;

		m_concentrations[ctr] = fC;

		ctr++;
	}

	return m_concentrations;
}

template <typename CAESReal>
typename TotalEquilibrium<CAESReal, true>::DDPack TotalEquilibrium<CAESReal, true>::distributionAnddDistdV(const CAESReal &v) const
{
	const size_t len = Ls.size();

	CAESReal X = 0.0;
	CAESReal dX = 0.0;
	const auto ts = Ts(v, X);	/* Recalculate Ts */
	const auto dts = dTsdV(v, dX);	/* Recalculate dTsdV */

	std::vector<CAESReal> distribution;
	std::vector<CAESReal> dDistdV;

	distribution.resize(len);
	dDistdV.resize(len);

	size_t ctr = 0;

	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = ts[idx];
		const CAESReal &dT = dts[idx];

		/* Distribution */
		const CAESReal fC = T / X;
		distribution[ctr] = fC;

		/* dDistdV */
		const CAESReal fD = (dT * X - T * dX) / (X * X);
		dDistdV[idx] = fD;
	}

	return { distribution, dDistdV };
}

template <typename CAESReal>
typename TotalEquilibrium<CAESReal, false>::DDPack TotalEquilibrium<CAESReal, false>::distributionAnddDistdV(const CAESReal &v)
{
	const size_t len = Ls.size();

	CAESReal X = 0.0;
	CAESReal dX = 0.0;
	Ts(v, X);	/* Recalculate Ts */
	dTsdV(v, dX);	/* Recalculate dTsdV */

	size_t ctr = 0;

	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = m_Ts[idx];
		const CAESReal &dT = m_dTsdV[idx];

		/* Distribution */
		const CAESReal fC = T / X;
		m_distribution[ctr] = fC;

		/* dDistdV */
		const CAESReal fD = (dT * X - T * dX) / (X * X);
		m_dDistdV[idx] = fD;
	}

	return { m_distribution, m_dDistdV };
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::distribution(const CAESReal &v) const
{
	CAESReal X = 0.0;
	const auto ts = Ts(v, X); /* Recalculate Ts */
	std::vector<CAESReal> distribution;

	distribution.resize(ts.size());
	size_t ctr = 0;
	for (const CAESReal &T : ts) {
		const CAESReal fC = T / X;

		distribution[ctr] = fC;

		ctr++;
	}

	return distribution;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::distribution(const CAESReal &v)
{
	CAESReal X = 0.0;
	Ts(v, X); /* Recalculate Ts */

	size_t ctr = 0;
	for (const CAESReal &T : m_Ts) {
		const CAESReal fC = T / X;

		m_distribution[ctr] = fC;

		ctr++;
	}

	return m_distribution;
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::dTsdV(const CAESReal &v, CAESReal &X) const
{
	const size_t len = Ls.size();
	X = 0.0;
	std::vector<CAESReal> dts{};

	assert(len == static_cast<size_t>(numHigh - numLow) + 1);

	dts.resize(len);
	dts[0] = 0.0;
	CAESReal vPow = 1.0;
	for (size_t idx = 1; idx < len; idx++) {
		const CAESReal &T = idx * Ls[idx] * vPow;

		vPow *= v;

		dts[idx] = T;
		X += T;
	}

	return dts;
}

template <typename CAESReal>
std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::dTsdV(const CAESReal &v, CAESReal &X)
{
	const size_t len = m_dTsdV.size();
	X = 0.0;

	assert(len == static_cast<size_t>(numHigh - numLow) + 1);

	m_dTsdV[0] = 0;
	CAESReal vPow = 1.0;
	for (size_t idx = 1; idx < len; idx++) {
		const CAESReal &T = idx * Ls[idx] * vPow;

		vPow *= v;

		m_dTsdV[idx] = T;
		X += T;
	}

	return m_dTsdV;
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::Ts(const CAESReal &v, CAESReal &X) const
{
	const size_t len = Ls.size();
	X = 0.0;
	std::vector<CAESReal> ts{};

	assert(len == static_cast<size_t>(numHigh - numLow + 1));

	ts.resize(len);
	CAESReal vPow = 1.0;
	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = Ls[idx] * vPow;

		vPow *= v;

		ts[idx] = T;
		X += T;
	}

	return ts;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::Ts(const CAESReal &v, CAESReal &X)
{
	const size_t len = m_Ts.size();
	X = 0.0;

	assert(len == static_cast<size_t>(numHigh - numLow + 1));

	CAESReal vPow = 1.0;
	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = Ls[idx] * vPow;

		vPow *= v;

		m_Ts[idx] = T;
		X += T;
	}

	return m_Ts;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_TOTALEQUILIBRIUM_HPP


