#include "totalequilibrium.h"

namespace ECHMET {
namespace CAES {

template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 const std::vector<double> &activityCoefficients,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len,
			 const int numLow)
{
	X = 1.0;
	ts[0] = 1.0;
	dX = 0.0;
	dts[0] = 0.0;

	const double activityLow = activityCoefficients[std::abs(numLow)];
	const double vMul = v * activityCoefficients[1];
	int num = numLow + 1;
	double dvPow = 1.0;
	double vPow = vMul;
	for (size_t idx = 1; idx < len; idx++) {
		const double actTerm = activityLow / activityCoefficients[std::abs(num++)];
		const double LA = Ls[idx] * actTerm;

		const double T = LA * vPow;
		const double dT = idx * LA * dvPow;

		dvPow = vPow;
		vPow *= vMul;

		ts[idx] = T;
		X += T;

		dts[idx] = dT;
		dX += dT;
	}
}

template <>
void calculateTsAnddTsdV<double, InstructionSet::AVX>
			(double *const ECHMET_RESTRICT_PTR ts,
			 double *const ECHMET_RESTRICT_PTR dts,
			 const double &v,
			 double &X, double &dX,
			 const std::vector<double> &Ls, const size_t len)
{
	X = 1.0;
	ts[0] = 1.0;
	dX = 0.0;
	dts[0] = 0.0;

	const double vMul = v;
	double dvPow = 1.0;
	double vPow = vMul;
	for (size_t idx = 1; idx < len; idx++) {
		const double LA = Ls[idx];

		const double T = LA * vPow;
		const double dT = idx * LA * dvPow;

		dvPow = vPow;
		vPow *= vMul;

		ts[idx] = T;
		X += T;

		dts[idx] = dT;
		dX += dT;
	}
}

} // namespace CAES
} // namespace ECHMET
