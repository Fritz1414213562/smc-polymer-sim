#ifndef POLYMER_MDMC_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define POLYMER_MDMC_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include "src/ff/HarmonicBondForceFieldGenerator.hpp"
#include "src/ff/HarmonicAngleForceFieldGenerator.hpp"
#include "src/ff/ExcludedVolumeForceFieldGenerator.hpp"
#include "src/ff/PolynomialRepulsiveForceFieldGenerator.hpp"
#include "src/ff/PullingForceFieldGenerator.hpp"

#include "src/Topology.hpp"
#include "src/util/Utility.hpp"

#include <OpenMM.h>

// -----------------------------------------------------------------------------
// read local force field

HarmonicBondForceFieldGenerator
read_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic);

HarmonicAngleForceFieldGenerator
read_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic);

// ----------------------------------------------------------------------------
// read global force field

std::vector<std::pair<std::size_t, std::size_t>>
read_ignore_molecule_and_particles_within(const toml::value& ignore_table, const Topology& topology);

std::vector<std::pair<std::string, std::string>>
read_ignore_group(const toml::value& ignore_table);

ExcludedVolumeForceFieldGenerator
read_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

PolynomialRepulsiveForceFieldGenerator
read_polynomial_repulsive_ff_generator(
	const toml::value& global_ff_data, const std::size_t system_size,
	const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
	const double temperature, const bool use_periodic);


// -----------------------------------------------------------------------------
// read external force field
PullingForceFieldGenerator
read_pulling_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic);

#endif // POLYMER_MDMC_READ_TOML_FORCE_FIELD_GENERATOR_HPP
