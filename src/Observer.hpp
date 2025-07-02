#ifndef POLYMER_MCMD_OBSERVER_HPP
#define POLYMER_MCMD_OBSERVER_HPP

#include "SystemGenerator.hpp"

#include <string>
#include <fstream>
#include <OpenMM.h>

class ObserverBase
{
	public:
		virtual ~ObserverBase() = default;
		virtual void initialize(const std::unique_ptr<OpenMM::System>&) = 0;
		virtual void output(const std::size_t, const OpenMM::Context&) = 0;
		virtual void finalize() const = 0;
		virtual std::string name() const = 0;
};

class DCDObserver final : public ObserverBase
{
	public:
		DCDObserver(const std::string& file_prefix, const std::size_t total_step,
			const std::size_t save_interval, const float delta_t, const bool use_periodic);

		void initialize(const std::unique_ptr<OpenMM::System>& system_ptr) override;

		void output(const std::size_t step, const OpenMM::Context& context) override;

		void finalize() const override { return; }

		std::string name() const override { return "DCDObserver"; }

	private:
		std::string        pos_filename_;
		std::string        vel_filename_;
		std::size_t        total_step_;
		float              delta_t_;
		std::size_t        save_interval_;
		bool               use_periodic_;
		std::vector<float> buffer_pos_x_;
		std::vector<float> buffer_pos_y_;
		std::vector<float> buffer_pos_z_;
		std::vector<float> buffer_vel_x_;
		std::vector<float> buffer_vel_y_;
		std::vector<float> buffer_vel_z_;

	private:
		void write_dcd_header(
			std::ofstream& ofs, const std::unique_ptr<OpenMM::System>& system_ptr) const;
		void write_unitcell(std::ofstream& ofs, const OpenMM::State& state) const;
		void write_dcd_frame(std::ofstream& ofs, const OpenMM::State& state);
		void write_dcd_velocity(std::ofstream& ofs, const OpenMM::State& state);
};


class EnergyObserver final : public ObserverBase
{
	public:
		EnergyObserver(const std::string& file_prefix, const SystemGenerator& system_gen);

		void initialize(const std::unique_ptr<OpenMM::System>&) override;

		void output(const std::size_t step, const OpenMM::Context& context) override;

		void finalize() const override { return ; }

		std::string name() const override { return "Energy Observer"; }

	private:

		void write_energy(std::ofstream& ofs, int frame_num,
						  const OpenMM::Context& context) const;

	private:
		std::string                        ene_filename_;
		std::map<std::string, std::size_t> ffname_groupid_map_;
};


//class MonteCarloResultObserver final : public ObserverBase
//{
//	public:
//		MonteCarloResultObserver(const std::string& file_prefix, const SystemGenerator& system_gen);
//		void initialize(const std::unique_ptr<OpenMM::System>&) override;
//		void output(const std::size_t step, const OpenMM::Context& context) override;
//		void finalize() const override { return; }
//		std::string name() const override { return "MonteCarloResultObserver"; }
//
//};

#endif // POLYMER_MCMD_OBSERVER_HPP
