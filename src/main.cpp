#include <OpenMM.h>
#include <toml.hpp>
#include "src/util/Utility.hpp"
#include "src/util/Macro.hpp"
#include "src/util/Logger.hpp"
#include "src/input/ReadInput.hpp"
#include "Simulator.hpp"
#include <string>
#include <cmath>

Simulator parse_command_line(int argc, char** argv)
{
	if (argc == 2)
	{
		const std::string file_suffix = Utility::get_file_suffix(std::string(argv[1]));
		if (file_suffix == ".toml") return read_input(std::string(argv[1]));
		else
			throw std::runtime_error(
				"The file format, " + file_suffix + " is not supported.");
	}
	log_fatal("Usage: {} <input.toml>", std::string(argv[0]));
}

void simulate(Simulator& simulator)
{
	const auto start = std::chrono::system_clock::now();
	log_info("Start!");

	simulator.run();

    const auto stop  = std::chrono::system_clock::now();
    const auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    if(total < 1000.0)
    {
        log_info("elapsed time : {} [msec]", total);
    }
    else if(total < 1000.0 * 60.0)
    {
        log_info("elapsed time : {:.1f} [sec]", total * 0.001);
    }
    else if(total < 1000.0 * 60.0 * 60.0)
    {
        log_info("elapsed time : {:.1f} [min]", total * 0.001 * 0.0167);
    }
    else if(total < 1000.0 * 60.0 * 60.0 * 24.0)
    {
        log_info("elapsed time : {:.1f} [hr]", total * 0.001 * 0.0167 * 0.0167);
    }
    else
    {
        log_info("elapsed time : {:.1f} [day]", total * 0.001 * 0.0167 * 0.0167 * 0.0417);
    }
}


int main(int argc, char** argv)
{
    // dump library information
    log_info("OpenMM Library Information");
    log_info("    version                  : {}", OpenMM::Platform::getOpenMMVersion());
    log_info("    GPU platform plugin path : {}", GPU_PLATFORM_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    // Load any shared libraries containing GPU implementations
    log_info("loading OpenMM plugins...");
    OpenMM::Platform::loadPluginsFromDirectory(
            GPU_PLATFORM_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    for(const auto& error : OpenMM::Platform::getPluginLoadFailures())
    {
        log_warn(error);
    }

    try {
        Simulator simulator(parse_command_line(argc, argv));
        simulate(simulator);
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("%s\n", e.what());
        return 1; // failure!
    }
}
