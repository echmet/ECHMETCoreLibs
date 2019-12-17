#include <echmetelems.h>
#include <iostream>

std::string isSIMDAvailable(const std::string &name, bool ECHMET::CPUSIMD::* dp, ECHMET::CPUSIMD &simd)
{
	if (simd.*dp)
		return name + ": Available";
	return name + ": Not available";
}

int main()
{
	auto cpuIdent = ECHMET::cpuIdentifier();
	auto simd = ECHMET::cpuSupportedSIMD();

	if (cpuIdent != nullptr)
		std::cout << "CPU Identifier: " << cpuIdent->c_str() << std::endl;

	std::cout << isSIMDAvailable("SSE2", &ECHMET::CPUSIMD::SSE2, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE3", &ECHMET::CPUSIMD::SSE3, simd) << std::endl;
	std::cout << isSIMDAvailable("SSSE3", &ECHMET::CPUSIMD::SSSE3, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE41", &ECHMET::CPUSIMD::SSE41, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE42", &ECHMET::CPUSIMD::SSE42, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX", &ECHMET::CPUSIMD::AVX, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX2", &ECHMET::CPUSIMD::AVX2, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX512", &ECHMET::CPUSIMD::AVX512, simd) << std::endl;
	std::cout << isSIMDAvailable("FMA3", &ECHMET::CPUSIMD::FMA3, simd) << std::endl;

	return 0;
}
	
