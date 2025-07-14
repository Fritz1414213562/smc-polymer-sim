#include "ReadForceFieldGenerator.hpp"

#include "src/util/Utility.hpp"
#include "src/util/Constants.hpp"

#include <toml.hpp>

// -----------------------------------------------------------------------------
// read local force field

// HarmonicBond input is like below
// [[forcefields.local]]
// interaction = "BondLength"
// potential   = "Harmonic"
// topology    = "bond"
// parameters  = [
//     {indices = [ 0, 1], v0 = 1.0, k = 2.0},
//     {indices = [ 3, 5], v0 = 3.0, k = 4.0},
//     {indices = [ 2, 6], v0 = 5.0, k = 6.0}
//     # ...
// ]
HarmonicBondForceFieldGenerator
read_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
		const double temperature, const bool use_periodic)
{
    Utility::check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(indices, offset);

        if(topology.size() <= indices.first)
        {
            throw std::runtime_error("[error] read_harmonic_bond_ff_generator : index "+std::to_string(indices.first)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }
        else if(topology.size() <= indices.second)
        {
            throw std::runtime_error("[error] read_harmonic_bond_ff_generator : index "+std::to_string(indices.second)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double v0 =
            Utility::find_parameter<double>(param, env, "v0");// * OpenMM::NmPerAngstrom; // nm
        // Toml input file assume the potential formula of HarmonicBond is
        // "k*(r - r0)^2", but OpenMM HarmonicBond is "1/2*k*(r - r0)^2".
        // So we needs double the interaction coefficient `k`.
        const double k =
            Utility::find_parameter<double>(param, env, "k")// * OpenMM::KJPerKcal * 2.0; // KJ/mol
				* 2.0
				* Constant::kB // J/K
				* Constant::Na // J/(K * mol)
				* 1e-3         // kJ / (K * mol)
				* temperature; // kJ / mol

            //OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * 2.0; // KJ/(mol nm^2)

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : Harmonic (" << indices_vec.size() << " found)" << std::endl;
    return HarmonicBondForceFieldGenerator(indices_vec, v0s, ks, use_periodic);
}

// HamonicAngle input is like below
// [[forcefields.local]]
// interaction   = "BondAngle"
// potential     = "Harmonic"
// topology      = "angle"
// parameters  = [
//     {indices = [ 0, 1, 2], k = 1.0, v0 = 2.0},
//     {indices = [ 3, 4, 5], k = 3.0, v0 = 4.0},
//     {indices = [ 2, 6, 7], k = 5.0, v0 = 6.0}
//     # ...
// ]
HarmonicAngleForceFieldGenerator
read_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
		const double temperature, const bool use_periodic)
{
    Utility::check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     v0s;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 3>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(indices, offset);

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_harmonic_angle_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }

        double v0 =
            Utility::find_parameter<double>(param, env, "v0"); // radian
        if(v0 < 0 || Constant::pi*2 < v0)
        {
            throw std::runtime_error("[error] read_harmonic_angle_ff_generator: "
                "v0 must be between 0 and 2pi");
        }
        else if(Constant::pi < v0)
        {
            v0 = 2.0 * Constant::pi - v0;
        }

        // Toml input file assume the potential formula of HarmonicAngle is
        // "k*(theta - theta0)^2", but OpenMM HarmonicAngle is
        // "1/2*k*(theta - theta0)^2". So we needs double the interaction
        // coefficient `k`.
        const double k =
            Utility::find_parameter<double>(param, env, "k")// * OpenMM::KJPerKcal * 2.0; // KJ/mol
				* 2.0
				* Constant::kB // J/K
				* Constant::Na // J/(K * mol)
				* 1e-3         // kJ / (K * mol)
				* temperature; // kJ / mol


        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    std::cerr << "    BondAngle     : Harmonic (" << indices_vec.size() << " found)" << std::endl;

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    return HarmonicAngleForceFieldGenerator(indices_vec, v0s, ks, use_periodic);
}

// GrosbergAngle input is like below
// [[forcefields.local]]
// interaction   = "BondAngle"
// potential     = "Grosberg"
// topology      = "angle"
// parameters  = [
//     {indices = [ 0, 1, 2], k = 1.0, v0 = 2.0},
//     {indices = [ 3, 4, 5], k = 3.0, v0 = 4.0},
//     {indices = [ 2, 6, 7], k = 5.0, v0 = 6.0}
//     # ...
// ]
GrosbergAngleForceFieldGenerator
read_grosberg_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
		const double temperature, const bool use_periodic)
{
    Utility::check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     v0s;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 3>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(indices, offset);

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_harmonic_angle_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }

        double v0 =
            Utility::find_parameter<double>(param, env, "v0"); // radian
        if(v0 < 0 || Constant::pi*2 < v0)
        {
            throw std::runtime_error("[error] read_harmonic_angle_ff_generator: "
                "v0 must be between 0 and 2pi");
        }
        else if(Constant::pi < v0)
        {
            v0 = 2.0 * Constant::pi - v0;
        }

        // Toml input file assume the potential formula of HarmonicAngle is
        // "k*(theta - theta0)^2", but OpenMM HarmonicAngle is
        // "1/2*k*(theta - theta0)^2". So we needs double the interaction
        // coefficient `k`.
        const double k =
            Utility::find_parameter<double>(param, env, "k")// * OpenMM::KJPerKcal * 2.0; // KJ/mol
				* Constant::kB // J/K
				* Constant::Na // J/(K * mol)
				* 1e-3         // kJ / (K * mol)
				* temperature; // kJ / mol


        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    std::cerr << "    BondAngle     : Grosberg (" << indices_vec.size() << " found)" << std::endl;

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    return GrosbergAngleForceFieldGenerator(indices_vec, v0s, ks, use_periodic);
}

// Segment Parallelization input is like below
// [[forcefields.local]]
// interaction = "BondLength"
// potential   = "SegmentParallelization"
// topology    = "bond"
// parameters  = [
//     {indices = [ 0, 1, 10, 11], v0 = 1.0, bond_k =  2.0, theta0 =  3.0, dihedral_k =  4.0},
//     {indices = [ 4, 5, 24, 25], v0 = 5.0, bond_k =  6.0, theta0 =  7.0, dihedral_k =  8.0},
//     {indices = [ 8, 9, 38, 39], v0 = 9.0, bond_k = 10.0, theta0 = 11.0, dihedral_k = 12.0}
//     # ...
// ]
SegmentParallelizationForceFieldGenerator
read_segment_parallelization_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
		const double temperature, const bool use_periodic)
{
    Utility::check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                             v0s;
    std::vector<double>                         theta0s;
    std::vector<double>                         bond_ks;
    std::vector<double>                     dihedral_ks;
	std::vector<double>                           phi0s;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(indices, offset);

		for (const auto& index : indices)
		{
			if(topology.size() <= index)
        	{
        	    throw std::runtime_error("[error] read_segment_parallelization_ff_generator : "
					"index " + std::to_string(index)
					+ " exceeds the system's largest index "
					+ std::to_string(topology.size()-1)+".");
        	}
        }

        const double v0 =
            Utility::find_parameter<double>(param, env, "v0");
		const double theta0 =
			Utility::find_parameter<double>(param, env, "theta0");
        // Toml input file assume the potential formula of HarmonicBond is
        // "k*(r - r0)^2", but OpenMM HarmonicBond is "1/2*k*(r - r0)^2".
        // So we needs double the interaction coefficient `k`.
        const double bond_k =
            Utility::find_parameter<double>(param, env, "bond_k")
				* Constant::kB // J/K
				* Constant::Na // J/(K * mol)
				* 1e-3         // kJ / (K * mol)
				* temperature; // kJ / mol
		const double dihedral_k =
			Utility::find_parameter<double>(param, env, "dihedral_k")
				* Constant::kB
				* Constant::Na
				* 1e-3
				* temperature;
		const double phi0 =
			Utility::find_parameter<double>(param, env, "phi0");

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        theta0s    .push_back(theta0);
        bond_ks    .push_back(bond_k);
        dihedral_ks.push_back(dihedral_k);
		phi0s      .push_back(phi0);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : SegmentParallelization ("
			  << indices_vec.size() << " found)" << std::endl;
    return SegmentParallelizationForceFieldGenerator(
		indices_vec, bond_ks, dihedral_ks, v0s, phi0s, theta0s, use_periodic);
}


// ----------------------------------------------------------------------------
// read global force field

std::vector<std::pair<std::size_t, std::size_t>>
read_ignore_molecule_and_particles_within(const toml::value& ignore_table, const Topology& topology)
{
    std::vector<std::pair<std::size_t, std::size_t>> ignore_list;
    bool                                             ignore_molecule_flag = false;

    if(ignore_table.contains("molecule"))
    {
        const std::string name = toml::find<std::string>(ignore_table, "molecule");
        if(name == "Intra" || name == "Self")
        {
            ignore_molecule_flag = true;
            std::vector<std::pair<std::size_t, std::size_t>> mol_ignore_list =
                topology.ignore_list_within_molecule();
            ignore_list.insert(ignore_list.end(),
                 mol_ignore_list.begin(), mol_ignore_list.end());
        }
        else if(name == "Others" || name == "Inter")
        {
            throw std::runtime_error(
                "[error] ignore molecule do not support \"Others\" or \"Inter\".");
        }
    }

    if(ignore_table.contains("particles_within"))
    {
        if(ignore_molecule_flag)
        {
            std::cerr << "\033[33m[warning]\033[m"
                      << " ignore molecule \"Intra\" or \"Self\" was defined,"
                         "so this ignore particle within bond will be ignored."
                      << std::endl;
        }

        const auto particle_within =
            toml::find<std::map<std::string, std::size_t>>(ignore_table, "particles_within");
        for(const auto& connection : particle_within)
        {
            std::vector<std::pair<std::size_t, std::size_t>> additional_list
                = topology.ignore_list_within_edge(connection.second, connection.first);
            ignore_list.insert(ignore_list.end(),
                    additional_list.begin(), additional_list.end());
        }
    }

    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    return ignore_list;
}

std::vector<std::pair<std::string, std::string>>
read_ignore_group(const toml::value& ignore_table)
{
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;

    if(ignore_table.contains("group"))
    {
        const auto& group = toml::find(ignore_table, "group");
        if(group.contains("inter"))
        {
            const auto& group_pairs = toml::find<toml::array>(group, "inter");
            for(const auto& group_pair : group_pairs)
            {

                const auto pair = toml::get<std::array<std::string, 2>>(group_pair);
                std::cerr << "        interaction between " << pair[0] << " and " << pair[1]
                           << " will be ignored" << std::endl;
                ignore_group_pairs.push_back({ pair[0], pair[1] });
            }
        }

        if(group.contains("intra"))
        {
            const auto&  groups = toml::find<toml::array>(group, "intra");
            for(const auto& group : groups)
            {
                const std::string group_str = toml::get<std::string>(group);
                std::cerr << "        ignore group intra in "
                          << group_str << " specified" << std::endl;

                ignore_group_pairs.push_back(std::make_pair(group_str, group_str));
            }
        }
    }

    return ignore_group_pairs;
}

// ExcludedVolume input is like below
// [[forcefield.global]]
// interaction                     = "Pair"
// potential                       = "ExcludedVolume"
// ignore.molecule                 = "Nothing"
// ignore.particles_within.bond    = 3
// ignore.particles_within.contact = 1
// epsilon                         = 0.6
// cutoff                          = 0.5
// parameters = [
// {index = 0, radius = 1.0},
// # ...
// ]
ExcludedVolumeForceFieldGenerator
read_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic)
{
    Utility::check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "env",
             "cutoff", "parameters", "epsilon"});

    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const double eps =
        toml::find<double>(global_ff_data, "epsilon") * OpenMM::KJPerKcal; // KJPermol
    const double cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_excluded_volume_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double      radius =
            Utility::find_parameter<double>(param, env, "radius");// * OpenMM::NmPerAngstrom; // nm
        radius_vec[index] = radius;
    }

    std::cerr << "    Global        : ExcludedVolume (" << params.size()
               << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return ExcludedVolumeForceFieldGenerator(
            eps, cutoff, radius_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

// ExcludedVolume input is like below
// [[forcefield.global]]
// interaction                     = "Pair"
// potential                       = "PolynomialRepulsive"
// ignore.molecule                 = "Nothing"
// ignore.particles_within.bond    = 1
// ignore.particles_within.contact = 1
// epsilon                         = 0.6
// cutoff                          = 1.0
// parameters = [
// {index = 0, radius = 1.0},
// # ...
// ]
PolynomialRepulsiveForceFieldGenerator
read_polynomial_repulsive_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const double temperature, const bool use_periodic)
{
    Utility::check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "env",
             "cutoff", "parameters", "epsilon"});

    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const double eps =
        toml::find<double>(global_ff_data, "epsilon")
		* Constant::kB // J/K
		* Constant::Na // J/(K * mol)
		* 1e-3         // kJ / (K * mol)
		* temperature; // kJ / mol

    const double cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_excluded_volume_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double      radius =
            Utility::find_parameter<double>(param, env, "radius");// * OpenMM::NmPerAngstrom; // nm
        radius_vec[index] = radius;
    }

    std::cerr << "    Global        : PolynomialRepulsive (" << params.size()
               << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return PolynomialRepulsiveForceFieldGenerator(
            eps, cutoff, radius_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

// -----------------------------------------------------------------------------
// read external force field

// PullingForce input is like below
// [[forcefields.external]]
// interaction = "PullingForce"
// parameters = [
//     {index = 0, force = [1.0, 0.0, 0.0]},
//     {index = 1, force = [0.0, 1.0, 0.0]},
//     # ...
// ]
PullingForceFieldGenerator
read_pulling_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic)
{
    using parameter_type = PullingForceFieldGenerator::parameter_type;

    Utility::check_keys_available(external_ff_data, {"interaction", "parameters", "env"});

    const auto& params = toml::find<toml::array>(external_ff_data, "parameters");
    const auto& env =
        external_ff_data.contains("env") ? external_ff_data.at("env") : toml::value{};

    std::vector<parameter_type> idx_force_vec;
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(index, offset);
        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_pulling_ff_generator : index " +
                std::to_string(index) + " exceeds the system's largest index " +
                std::to_string(topology.size()-1) + ".");
        }
        const std::array<double, 3> force_kcal_vec =
            Utility::find_parameter<std::array<double, 3>>(param, env, "force"); // [kcal/(mol Å)]
        std::array<double, 3> force_kj_vec;
        for(std::size_t idx = 0; idx < 3; idx++)
        {
            force_kj_vec[idx] =
                force_kcal_vec[idx] * OpenMM::KJPerKcal;// * OpenMM::AngstromsPerNm;
        }
        idx_force_vec.push_back({index, force_kj_vec});
    }

    std::cerr << "    External      : Pulling (" << params.size()
              << " found)" << std::endl;

    return PullingForceFieldGenerator(idx_force_vec, use_periodic);
}

// PositionRestraint input is like below
// [[forcefields.external]]
// interaction = "PositionRestraint"
// potential   = "Harmonic"
// parameters = [
//     {index = 0, position = [0.0, 0.0, 0.0], k = 0.1, v0 = 10.0},
//     # ...
// ]
PositionRestraintForceFieldGenerator
read_position_restraint_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic)
{
    Utility::check_keys_available(external_ff_data,
            {"interaction", "potential", "parameters", "env"});

    const auto& params = toml::find<toml::array>(external_ff_data, "parameters");
    const auto& env =
        external_ff_data.contains("env") ? external_ff_data.at("env") : toml::value{};

    std::vector<std::size_t>           indices;
    std::vector<std::array<double, 3>> positions;
    std::vector<double>                ks;
    std::vector<double>                v0s;
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        Utility::add_offset(index, offset);
        if(topology.size() <= index)
        {
            throw std::runtime_error(
                "[error] read_toml_position_restraint_ff_generator : index " +
                std::to_string(index) + " exceeds the system's largest index " +
                std::to_string(topology.size()-1) + ".");
        }
        const std::array<double, 3> position =
            Utility::find_parameter<std::array<double, 3>>(param, env, "position"); // [nm]
        const double k =
            Utility::find_parameter<double>(param, env, "k" );
        const double v0 =
            Utility::find_parameter<double>(param, env, "v0"); // [Nm]

        indices  .push_back(index);
        positions.push_back(position);
        ks       .push_back(k);
        v0s      .push_back(v0);
    }

    std::cerr << "    External      : PositionRestraint (" << params.size()
              << " found)" << std::endl;

    return PositionRestraintForceFieldGenerator(indices, positions, ks, v0s, use_periodic);
}

