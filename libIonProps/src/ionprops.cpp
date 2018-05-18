#include "ionprops_p.h"

namespace ECHMET {
namespace IonProps {

ECHMETReal ECHMET_CC calculateConductivity(const ComputationContext *ctx, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculateConductivityInternal<ECHMETReal>(ctx, calcProps);
}

RetCode ECHMET_CC calculateEffectiveMobilities(const ComputationContext *ctx, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculateEffectiveMobilitiesInternal<ECHMETReal>(ctx, analyticalConcentrations, calcProps);
}

ECHMETReal ECHMET_CC calculatepH(const ComputationContext *ctx, const NonidealityCorrections corrections, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculatepHInternal<ECHMETReal>(ctx, corrections, calcProps);
}

ECHMETReal ECHMET_CC calculatepH_direct(const ECHMETReal &cH, const ECHMETReal &ionicStrength) noexcept
{
	return calculatepH_directInternal<ECHMETReal>(cH, ionicStrength);
}

RetCode ECHMET_CC correctMobilities(const ComputationContext *ctx, const NonidealityCorrections corrections, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	return correctMobilitiesInternal<ECHMETReal>(ctx, corrections, analyticalConcentrations, calcProps);
}

ComputationContext * ECHMET_CC makeComputationContext(const SysComp::ChemicalSystem &chemSystem) noexcept
{
	return new (std::nothrow) ComputationContextImpl<ECHMETReal>(chemSystem);
}

ComputationContext::~ComputationContext() noexcept {}

} // namespace IonProps

} // namespace ECHMET
