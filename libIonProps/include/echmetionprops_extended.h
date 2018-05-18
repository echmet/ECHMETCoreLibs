#ifndef ECHMET_IONPROPS_IONPROPS_EXTENDED_H
#define ECHMET_IONPROPS_IONPROPS_EXTENDED_H

#include "../src/ionprops_p.h"

namespace ECHMET {
namespace IonProps {

template <bool B, typename IPReal>
struct ContextMaker {
	//static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const std::vector<IPReal> &ionicConcentrations);
};

template <typename IPReal>
struct ContextMaker<true, IPReal> {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations)
	{
		(void)ionicConcentrations;
		return new (std::nothrow) ComputationContextImpl<IPReal>{chemSystem};
	}
};

template <typename IPReal>
struct ContextMaker<false, IPReal> {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations)
	{
		return new (std::nothrow) ComputationContextImpl<IPReal>{ionicConcentrations, chemSystem};
	}
};

template <typename IPReal>
ComputationContext * makeComputationContextExtended(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations)
{
	return ContextMaker<std::is_same<IPReal, ECHMETReal>::value, IPReal>::make(chemSystem, ionicConcentrations);
}

} // namespace IonProps
} // namespace ECHMET


#endif // ECHMET_IONPROPS_IONPROPS_EXTENDED_H
