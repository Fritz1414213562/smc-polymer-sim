#include "IntegratorGenerator.hpp"
#include <fmt/core.h>


template<typename IntegratorType>
std::unique_ptr<OpenMM::Integrator> IntegratorGenerator::generate() const
{
	auto integrator = std::make_unique<IntegratorType>(
		this->temperature_,
		this->friction_coeff_,
		this->delta_t);
	integrator->setRandomNumberSeed(this->seed_);
	return integrator;
}

std::string IntegratorGenerator::dump_info() const
{
	std::string info;
	info += fmt::format("    temperature           : {:10.2f} K\n",     this->temperature_);
	info += fmt::format("    friction_gamma        : {:10.3f} ps^-1\n", this->friction_coeff_);
	info += fmt::format("    delta_t (or tolerance): {:10.3f} ps\n", this->delta_t_);

	std::string seed_info;
	if (this->seed_ == 0)
	{
		seed_info = "not specified or 0. random seed will be chosen.";
	}
	else
	{
		seed_info = fmt::format("{:10}", this->seed_);
	}
	info += fmt::format("    seed                  : {}\n", seed_info);
	return info;
}
