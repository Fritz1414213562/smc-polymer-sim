#include "ReadTOMLInput.hpp"

#include "src/util/Utility.hpp"
#include "src/util/Logger.hpp"

#include "src/Simulator.hpp"
#include "src/SystemGenerator.hpp"
#include "src/IntegratorGenerator.hpp"
#include "src/Topology.hpp"
#include "ReadForceFieldGenerator.hpp"
#include "src/ff/MonteCarloAnisotropicBarostatGenerator.hpp"

#include <OpenMM.h>
#include <toml.hpp>

#include <memory>
#include <iostream>

SystemGenerator read_system(const toml::value& data)
{
    // read systems table
    const auto& systems   = toml::find(data, "systems");

    // read particles info
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    Topology        topology(system_size);
    std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    std::vector<double>                     mass_vec(system_size);
    for(std::size_t idx=0; idx<system_size; ++idx)
    {
        // set mass
        const auto& p = particles.at(idx);
        mass_vec[idx] = toml::get<double>(Utility::find_either(p, "m", "mass")); // amu

        if(p.contains("group"))
        {
            const std::string group_name = toml::find<std::string>(p, "group");
            group_vec[idx] = group_name;
        }
    }
    SystemGenerator system_gen(mass_vec);

    // read boundary condition
    const auto&       simulator_table = toml::find(data, "simulator");
    const std::string boundary_type   = toml::find<std::string>(simulator_table, "boundary_type");
    bool use_periodic = false;
    if(boundary_type == "Periodic" || boundary_type == "PeriodicCuboid")
    {
        use_periodic = true;
        const auto& boundary_shape = toml::find(systems[0], "boundary_shape");
        const auto lower_bound = toml::find<std::array<double, 3>>(boundary_shape, "lower");
        const auto upper_bound = toml::find<std::array<double, 3>>(boundary_shape, "upper");
        if(lower_bound[0] != 0.0 || lower_bound[1] != 0.0 || lower_bound[2] != 0.0)
        {
            std::cerr << "\033[33m[warning]\033[m"
                      << " Lower bound of periodic boundary box is not (0.0, 0.0, 0.0). "
                         "OpenAICG2+ only support the case lower bound is origin, "
                         "so the simulation box will be moved to satisfy this condition. "
                         "In this case, specified lower bound is ("
                      << std::fixed << std::setprecision(2)
                      << std::setw(7) << lower_bound[0] << ", "
                      << std::setw(7) << lower_bound[1] << ", "
                      << std::setw(7) << lower_bound[2] << ")"
                      << std::endl;
        }

        const double xlength = (upper_bound[0] - lower_bound[0]) * OpenMM::NmPerAngstrom;
        const double ylength = (upper_bound[1] - lower_bound[1]) * OpenMM::NmPerAngstrom;
        const double zlength = (upper_bound[2] - lower_bound[2]) * OpenMM::NmPerAngstrom;
        system_gen.set_pbc(xlength, ylength, zlength);
    }

    // read ensemble condition
    const auto& attr = toml::find(systems[0], "attributes");
    std::cerr << "reading ensemble condition..." << std::endl;
    if(attr.contains("ensemble"))
    {
        const auto& ensemble = toml::find(attr, "ensemble");
        const auto& type     = toml::find<std::string>(ensemble, "type");
        if(type == "NVT")
        {
            std::cerr << "    ensemble type is NVT" << std::endl;
        }
        else if(type == "NPT")
        {
            std::cerr << "    ensemble type is NPT with anisotropic barostat" << std::endl;
            const auto& default_pressure =
                toml::find<std::array<double, 3>>(ensemble, "pressure");
            const auto& scale_axis =
                toml::find<std::array<bool, 3>>(ensemble, "scale_axis");
            const auto frequency = toml::find_or<std::size_t>(ensemble, "frequency", 25);

            if(!attr.contains("temperature"))
            {
                throw std::runtime_error(
                        "[error] attributes table must contains temperature in NPT ensemble case.");
            }
            const double temperature = toml::find<double>(attr, "temperature");

            auto barostat_gen =
                MonteCarloAnisotropicBarostatGenerator(
                        scale_axis, temperature, default_pressure, frequency);
            system_gen.set_barostat(
                    std::make_unique<MonteCarloAnisotropicBarostatGenerator>(barostat_gen));
        }
        else
        {
            throw std::runtime_error(
                    "[error] unknown ensemble type "+type+" was specified.");
        }
    }

    std::cerr << "reading forcefield tables..." << std::endl;
    // read forcefields info
    const auto ff = toml::find(data, "forcefields").at(0);
    if(ff.contains("local"))
    {
        const auto& locals = toml::find(ff, "local").as_array();

        for(const auto& local_ff : locals)
        {
            const std::string interaction = toml::find<std::string>(local_ff, "interaction");
            const std::string potential   = toml::find<std::string>(local_ff, "potential");
            if(interaction == "BondLength" && potential == "Harmonic")
            {
                HarmonicBondForceFieldGenerator ff_gen =
                    read_harmonic_bond_ff_generator(
                            local_ff, topology, use_periodic);
                system_gen.add_ff_generator(
                        std::make_unique<HarmonicBondForceFieldGenerator>(ff_gen));
            }
            else if(interaction == "BondAngle" && potential == "Harmonic")
            {
                HarmonicAngleForceFieldGenerator ff_gen =
                    read_harmonic_angle_ff_generator(local_ff, topology, use_periodic);
                system_gen.add_ff_generator(
                        std::make_unique<HarmonicAngleForceFieldGenerator>(ff_gen));
            }
        }
    }
    topology.make_molecule("bond");

    if(ff.contains("global"))
    {
        const auto& globals = toml::find(ff, "global").as_array();

        for(const auto& global_ff : globals)
        {
            const std::string interaction =
                toml::find<std::string>(global_ff, "interaction");
            const std::string potential =
                toml::find<std::string>(global_ff, "potential");

            // Pair interaction case
            if(potential == "ExcludedVolume")
            {
                ExcludedVolumeForceFieldGenerator ff_gen =
                    read_excluded_volume_ff_generator(
                        global_ff, system_size, topology, group_vec, use_periodic);
                system_gen.add_ff_generator(
                        std::make_unique<ExcludedVolumeForceFieldGenerator>(ff_gen));
            }
			else if (potential == "PolynomialRepulsive")
			{
				if(!attr.contains("temperature"))
            	{
            	    throw std::runtime_error(
            	        "[error] attributes table must contains temperature for PolynomialRepulsive");
            	}
            	const double temperature = toml::find<double>(attr, "temperature");
				PolynomialRepulsiveForceFieldGenerator ff_gen =
					read_polynomial_repulsive_ff_generator(
						global_ff, system_size, topology, group_vec, temperature, use_periodic);
				system_gen.add_ff_generator(
						std::make_unique<PolynomialRepulsiveForceFieldGenerator>(ff_gen));
			}
        }
    }

    if(ff.contains("external"))
    {
        const auto& externals = toml::find(ff, "external").as_array();

        for(const auto& external_ff : externals)
        {
            const std::string interaction =
                toml::find<std::string>(external_ff, "interaction");
            if(interaction == "PullingForce")
            {
                PullingForceFieldGenerator ff_gen =
                    read_pulling_ff_generator(external_ff, topology, use_periodic);
                system_gen.add_ff_generator(
                    std::make_unique<PullingForceFieldGenerator>(ff_gen));
            }
        }
    }

    return system_gen;
}

std::unique_ptr<IntegratorGeneratorBase>
read_integrator_gen(const toml::value& root)
{
    // setup OpenMM integrator

    const auto& systems     = toml::find(root, "systems");
    const auto& attr        = toml::find(systems[0], "attributes");
    const auto  temperature = toml::find<double>(attr, "temperature");

    const auto& simulator   = toml::find(root, "simulator");
    const auto  delta_t     =
        toml::find<double>(simulator, "delta_t") * Constant::cafetime; // [ps]

    const auto  seed = toml::find_or<int>(simulator, "seed", 0);

    // In OpenMM, we cannot use different friction coefficiet, gamma, for
    // different molecules. However, cafemol use different gamma depends on the
    // mass of each particle, and the product of mass and gamma is constant in
    // there. So in this implementation, we difine default value of gamma to 0.2
    // ps^-1, correspond to approximatry 0.01 in cafemol friction coefficient.
    // We need to implement new LangevinIntegrator which can use different gamma
    // for different particles.

    const auto  gamma_t =
        toml::find_or<double>(simulator, "gamma_t", 0.01); // [1/cafetime]
    const double friction_coeff = gamma_t / Constant::cafetime; // [1/ps]

    return std::make_unique<LangevinIntegratorGenerator>(
            temperature, friction_coeff, delta_t, seed);
}

std::vector<OpenMM::Vec3> read_initial_conf(const toml::value& data)
{
    const auto& systems     = toml::find(data, "systems");
    const auto& particles   = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initPosInNm(system_size);
    for (std::size_t i=0; i<system_size; ++i)
    {
        // set position
        const auto& p = particles.at(i);
        std::array<double, 3> vec = {0.0, 0.0, 0.0};
        vec = toml::get<std::array<double, 3>>(Utility::find_either(p, "pos", "position"));
        initPosInNm[i] = OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                      vec[1]*OpenMM::NmPerAngstrom,
                                      vec[2]*OpenMM::NmPerAngstrom); // location, nm
    }

    return initPosInNm;
}

std::vector<OpenMM::Vec3> read_initial_vel(const toml::value& data)
{
    const auto& systems   = toml::find(data, "systems");
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initVelInNmPs(system_size); // [nm/ps]
    if(particles.at(1).contains("vel") || particles.at(1).contains("velocity"))
    {
        for(std::size_t i=0; i<system_size; ++i)
        {
            // set velocity
            const auto& p = particles.at(i);
            if(p.contains("vel") || p.contains("velocity"))
            {
                const auto vec =
                    toml::get<std::array<double, 3>>(
                            Utility::find_either(p, "vel", "velocity"));
                initVelInNmPs[i] =
                    OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                 vec[1]*OpenMM::NmPerAngstrom,
                                 vec[2]*OpenMM::NmPerAngstrom);
            }
            else
            {
                throw std::runtime_error(
                        "[error] All particle need to have velocity,"
                        " because first particle in system has velocity.");
            }
        }
    }

    return initVelInNmPs;
}

Simulator read_input(const std::string& file_name)
{
    std::size_t file_path_len  = file_name.rfind("/")+1;
    if(file_path_len == std::string::npos)
    {
        file_path_len = 0;
    }

    const std::string file_suffix = Utility::get_file_suffix(file_name);
    if(file_suffix != ".toml")
    {
            throw std::runtime_error(
                    "[error] File suffix is `" + file_suffix + "`."
                    " Toml mode needs `.toml` file for toml.");
    }

    // read toml toml file
    std::cerr << "parsing " << file_name << "..." << std::endl;
    auto data = toml::parse(toml_file_name);
    expand_include(data);

    // read files table
    const auto&        files         = toml::find(data, "files");
    const auto&        output        = toml::find(files, "output");
    const std::string& output_prefix = toml::find<std::string>(output, "prefix");
    const std::string& output_format = toml::find<std::string>(output, "format");
    const bool         dump_progress_bar =
        toml::find_or<bool>(output, "progress_bar", true);
    std::string output_path   = toml::find<std::string>(output, "path");
    if(output_path.empty())
    {
        output_path = "./";
    }
    if(output_path.back() != '/')
    {
        output_path += '/';
    }

    // read simulator table
    const auto&        simulator_table = toml::find(data, "simulator");
    const std::string& boundary_type   = toml::find<std::string>(simulator_table, "boundary_type");
    const std::size_t  total_step      = toml::find<std::size_t>(simulator_table, "total_step");
    const std::size_t  save_step       = toml::find<std::size_t>(simulator_table, "save_step");
    const double       delta_t         = toml::find<double>(simulator_table, "delta_t");
    const bool         energy_minimization =
        toml::find_or<bool>(simulator_table, "energy_minimization", false);

    // read system table
    const auto& systems     = toml::find(data, "systems");
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_initial_conf(data));
    SystemGenerator system_gen = read_system(data);

    std::unique_ptr<IntegratorGeneratorBase> integrator_gen_ptr =
        read_integrator_gen(data);

    // construct observers
    const bool use_periodic =
        (boundary_type == "Periodic" || boundary_type == "PeriodicCuboid");
    std::vector<std::unique_ptr<ObserverBase>> observers;
    if(output_format == "dcd")
    {
        observers.push_back(
                std::make_unique<DCDObserver>(
                    output_path+output_prefix, total_step, save_step, delta_t, use_periodic));
    }
    else
    {
        throw std::runtime_error(
                "[error] output file format `" + output_format + "` is not supported.");
    }
    observers.push_back(std::make_unique<EnergyObserver>(output_path+output_prefix, system_gen));

    // read platform
    const auto platform_table = toml::find_or(data, "platform", toml::value(toml::table{}));
    const auto platform_name  = toml::find_or<std::string>(platform_table, "name", std::string("CUDA"));
    std::cerr << "    OpenMM platform : " << platform_name << std::endl;

    // platform properties corresponds to OpenMM platform-specific properties.
    // The detail informations are in
    // http://docs.openmm.org/latest/userguide/library/04_platform_specifics.html
    const auto platform_properties =
        toml::find_or<std::map<std::string, std::string>>(platform_table, "properties", {/*no properties*/});

    // check if the platform is available
    bool platform_found = false;
    for(int i=0; i<OpenMM::Platform::getNumPlatforms(); ++i)
    {
        if(OpenMM::Platform::getPlatform(i).getName() == platform_name)
        {
            platform_found = true;
            break;
        }
    }
    if(!platform_found)
    {
        if(platform_table.contains("name"))
        {
            throw std::runtime_error(toml::format_error("[error] platform \"" +
                platform_name + "\" not found. You need to set the correct OpenMM "
                "plugins directory path to the CMake option -DOPENMM_PLUGIN_DIR.",
                platform_table.at("name"), "defined here"));
        }
        else
        {
            throw std::runtime_error("[error] platform \"" +
                platform_name + "\" not found. You need to set the correct OpenMM "
                "plugins directory path to the CMake option -DOPENMM_PLUGIN_DIR.");
        }
    }

    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName(platform_name);
    Simulator simulator =
        Simulator(system_gen, std::move(integrator_gen_ptr),
                  platform, platform_properties,
                  initial_position_in_nm, total_step, save_step,
                  observers, energy_minimization, dump_progress_bar);

    const auto& particles = toml::find<toml::array>(systems[0], "particles");
    if(particles.at(0).contains("vel") || particles.at(0).contains("velocity"))
    {
        const std::vector<OpenMM::Vec3>
            initial_vel_in_nmps(read_initial_vel(data));
        simulator.set_velocity(initial_vel_in_nmps);
    }
    else
    {
        const auto& sys = systems.at(0);
        const auto vel_seed = toml::find_or<std::int64_t>(sys, "vel_seed", 0);
        if(vel_seed == 0)
        {
            // since OpenMM does not provide system-wide RNG, we need to specify
            // seeds for each specific target.
            log_info("System does not provide initial velocity nor non-zero `vel_seed`. "
                     "initial velocity becomes completely random and cannot be reproduced.");
        }
        simulator.set_velocity(vel_seed);
    }

    return simulator;
}

