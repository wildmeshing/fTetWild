#pragma once

#include <floattetwild/Mesh.hpp>
#include <floattetwild/get_mem.h>
#include <mutex>

namespace floatTetWild {

    class StateInfo {
    public:
        static constexpr int tetgen_id = -1;

        static constexpr int init_id = 0;
        static constexpr int preprocessing_id = 1;
        static constexpr int tetrahedralization_id = 2;
        static constexpr int cutting_id = 3;
        static constexpr int optimization_id = 4;
        static constexpr int wn_id = 5;

		static constexpr int splitting_id = 6;
		static constexpr int collapsing_id = 7;
        static constexpr int swapping_id = 8;
        static constexpr int smoothing_id = 9;

        int id = 0;
        int v_num = 0;
        int t_num = 0;
        double time = 0;
        Scalar avg_energy = 0;
        Scalar max_energy = 0;
        int cnt_fail_inserted_face = -1;

        StateInfo(int _id, double _time, int _v_num, int _t_num, Scalar _max_energy, Scalar _avg_energy,
				  int _cnt_fail_inserted_face = -1):
        id(_id), time(_time), v_num(_v_num), t_num(_t_num), max_energy(_max_energy), avg_energy(_avg_energy),
		cnt_fail_inserted_face(_cnt_fail_inserted_face){}
    };

	class Statistics {
	public:
	    static Statistics & instance() {
	        Statistics inst;
	        return inst;
	    }

		template<class... Args>
	    void record(Args && ... args) {
	    	std::lock_guard<std::mutex> guard(mutex_);
	    	states_.emplace_back(std::forward<Args>(args)...);
	    }

		void print_sum() const {
			double time = 0;
			for (auto &s : states_) {
				if (s.id != StateInfo::splitting_id && s.id != StateInfo::collapsing_id
					&& s.id != StateInfo::swapping_id && s.id != StateInfo::smoothing_id)
					time += s.time;
			}
			cout << -1 << ", " << time << ", " << states_.back().v_num << ", " << states_.back().t_num << ", "
				 << states_.back().cnt_fail_inserted_face << endl;
		}

		friend std::ostream &operator<<(std::ostream &stream, const Statistics &stats) {
			double time = 0;
			int cnt_uninserted = 0;
			for (auto &s : stats.states_) {
				stream << s.id << ", " << s.time << ", " << s.v_num << ", " << s.t_num << ", "
					   << s.max_energy << ", " << s.avg_energy << ", " << s.cnt_fail_inserted_face << ", -1" << endl;
				if (s.cnt_fail_inserted_face >= 0)
					cnt_uninserted = s.cnt_fail_inserted_face;
				if (s.id < 6)
					time += s.time;
			}
			stream << -1 << ", " << time << ", " << stats.states_.back().v_num << ", "
				   << stats.states_.back().t_num << ", "
				   << stats.states_.back().max_energy << ", " << stats.states_.back().avg_energy << ", "
				   << cnt_uninserted << ", " << get_peak_mem() << endl;

			return stream;
		}

	private:
		Statistics() = default;

	private:
		std::mutex mutex_;
		std::vector<StateInfo> states_;
	};

	// Retrieve current stats instance
	inline Statistics & stats() {
		return Statistics::instance();
	}
}
