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
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations, const ComputationContext::Options options)
	{
		(void)ionicConcentrations;
		return new (std::nothrow) ComputationContextImpl<IPReal>{chemSystem, options};
	}
};

template <typename IPReal>
struct ContextMaker<false, IPReal> {
	static ComputationContext * make(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations, const ComputationContext::Options options)
	{
		return new (std::nothrow) ComputationContextImpl<IPReal>{ionicConcentrations, chemSystem, options};
	}
};

template <typename IPReal>
ComputationContext * makeComputationContextExtended(const SysComp::ChemicalSystem &chemSystem, const std::vector<IPReal> &ionicConcentrations, const ComputationContext::Options options)
{
	return ContextMaker<std::is_same<IPReal, ECHMETReal>::value, IPReal>::make(chemSystem, ionicConcentrations, options);
}

} // namespace IonProps
} // namespace ECHMET


#endif // ECHMET_IONPROPS_IONPROPS_EXTENDED_H
