#ifndef POLYMER_MDMC_READ_TOML_INPUT_HPP
#define POLYMER_MDMC_READ_TOML_INPUT_HPP

#include "src/Simulator.hpp"
#include "src/SystemGenerator.hpp"
#include "src/IntegratorGenerator.hpp"

#include <toml_fwd.hpp>
#include <OpenMM.h>

#include <memory>

SystemGenerator read_system(const toml::value& data);

std::unique_ptr<IntegratorGeneratorBase>
read_integrator_gen(const toml::value& root);

std::vector<OpenMM::Vec3> read_initial_conf(const toml::value& data);
std::vector<OpenMM::Vec3> read_initial_vel(const toml::value& data);

Simulator read_input(const std::string& file_name);

#endif // POLYMER_MDMC_READ_TOML_INPUT_HPP
