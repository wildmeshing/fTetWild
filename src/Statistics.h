#pragma once

#include <floattetwild/Mesh.hpp>
#include <floattetwild/get_mem.h>

namespace floatTetWild {
    class StateInfo{
    public:
        static const int tetgen_id = -1;

        static const int init_id = 0;
        static const int preprocessing_id = 1;
        static const int tetrahedralization_id = 2;
        static const int cutting_id = 3;
        static const int optimization_id = 4;
        static const int wn_id = 5;

		static const int splitting_id = 6;
		static const int collapsing_id = 7;
        static const int swapping_id = 8;
        static const int smoothing_id = 9;

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
		inline static Statistics &stats() {
			static Statistics stats;
			return stats;
		}

		std::vector<StateInfo> states;

		void print_sum() {
			double time = 0;
			for (auto &s:states) {
				if (s.id != StateInfo::splitting_id && s.id != StateInfo::collapsing_id
					&& s.id != StateInfo::swapping_id && s.id != StateInfo::smoothing_id)
					time += s.time;
			}
			cout << -1 << ", " << time << ", " << states.back().v_num << ", " << states.back().t_num << ", "
				 << states.back().cnt_fail_inserted_face << endl;
		}

		friend std::ostream &operator<<(std::ostream &stream, const Statistics &stats) {
			double time = 0;
			int cnt_uninserted = 0;
			for (auto &s:stats.states) {
				stream << s.id << ", " << s.time << ", " << s.v_num << ", " << s.t_num << ", "
					   << s.max_energy << ", " << s.avg_energy << ", " << s.cnt_fail_inserted_face << ", -1" << endl;
				if (s.cnt_fail_inserted_face >= 0)
					cnt_uninserted = s.cnt_fail_inserted_face;
				if (s.id < 6)
					time += s.time;
			}
			stream << -1 << ", " << time << ", " << stats.states.back().v_num << ", "
				   << stats.states.back().t_num << ", "
				   << stats.states.back().max_energy << ", " << stats.states.back().avg_energy << ", "
				   << cnt_uninserted << ", " << get_peak_mem() << endl;

			return stream;
		}

	private:
		Statistics() {}
	};
}
