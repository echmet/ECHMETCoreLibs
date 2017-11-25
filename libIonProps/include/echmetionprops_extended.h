#ifndef ECHMET_IONPROPS_IONPROPS_EXTENDED_H
#define ECHMET_IONPROPS_IONPROPS_EXTENDED_H

#include "../src/ionprops_p.h"

namespace ECHMET {
namespace IonProps {

template <bool B, typename IPReal>
struct ContextMaker {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const std::vector<IPReal> &ionicConcentrations);
};

template <typename IPReal>
struct ContextMaker<true, IPReal> {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const std::vector<IPReal> &ionicConcentrations)
	{
		(void)ionicConcentrations;
		return new (std::nothrow) ComputationContextImpl<IPReal>{chemSystem, analyticalConcentrations, calcProps};
	}
};

template <typename IPReal>
struct ContextMaker<false, IPReal> {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const std::vector<IPReal> &ionicConcentrations)
	{
		return new (std::nothrow) ComputationContextImpl<IPReal>{ionicConcentrations, chemSystem, analyticalConcentrations, calcProps};
	}
};

template <typename IPReal>
ComputationContext * makeComputationContextExtended(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps, const std::vector<IPReal> &ionicConcentrations)
{
	return ContextMaker<std::is_same<IPReal, ECHMETReal>::value, IPReal>::make(chemSystem, analyticalConcentrations, calcProps, ionicConcentrations);
}

} // namespace IonProps
} // namespace ECHMET


#endif // ECHMET_IONPROPS_IONPROPS_EXTENDED_H
