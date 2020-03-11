// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/get_mem.h>
#include <floattetwild/Mesh.hpp>
#include <mutex>

namespace floatTetWild {

    class StateInfo {
    public:
        enum OpId {
            tetgen_id = -1,

            init_id = 0,
            preprocessing_id = 1,
            tetrahedralization_id = 2,
            cutting_id = 3,
            optimization_id = 4,
            wn_id = 5,

            splitting_id = 6,
            collapsing_id = 7,
            swapping_id = 8,
            smoothing_id = 9
        };

        int id = 0;
        int v_num = 0;
        int t_num = 0;
        double time = 0;
        Scalar avg_energy = 0;
        Scalar max_energy = 0;
        int cnt_fail_inserted_face = -1;

        StateInfo(int _id,
                  double _time,
                  int _v_num,
                  int _t_num,
                  Scalar _max_energy,
                  Scalar _avg_energy,
                  int _cnt_fail_inserted_face = -1)
                : id(_id), time(_time), v_num(_v_num), t_num(_t_num), max_energy(_max_energy), avg_energy(_avg_energy),
                  cnt_fail_inserted_face(_cnt_fail_inserted_face) {}
    };

    class Statistics {
    public:
        static Statistics &instance() {
            static Statistics inst;
            return inst;
        }

        template<class... Args>
        void record(Args &&... args) {
            std::lock_guard<std::mutex> guard(mutex_);
            states_.emplace_back(std::forward<Args>(args)...);
        }

        void print_sum() const {
            double time = 0;
            for (auto &s : states_) {
                if (s.id != StateInfo::splitting_id && s.id != StateInfo::collapsing_id &&
                    s.id != StateInfo::swapping_id && s.id != StateInfo::smoothing_id)
                    time += s.time;
            }
            cout << -1 << ", " << time << ", " << states_.back().v_num << ", " << states_.back().t_num
                 << ", " << states_.back().cnt_fail_inserted_face << endl;
        }

        friend std::ostream &operator<<(std::ostream &stream, const Statistics &stats) {
            double time = 0;
            int cnt_uninserted = 0;
            for (auto &s : stats.states_) {
                stream << s.id << ", " << s.time << ", " << s.v_num << ", " << s.t_num << ", "
                       << s.max_energy << ", " << s.avg_energy << ", " << s.cnt_fail_inserted_face
                       << ", -1" << endl;
                if (s.cnt_fail_inserted_face >= 0)
                    cnt_uninserted = s.cnt_fail_inserted_face;
                if (s.id < 6)
                    time += s.time;
            }
            stream << -1 << ", " << time << ", " << stats.states_.back().v_num << ", "
                   << stats.states_.back().t_num << ", " << stats.states_.back().max_energy << ", "
                   << stats.states_.back().avg_energy << ", " << cnt_uninserted << ", "
                   << get_peak_mem() << endl;

            return stream;
        }

    private:
        Statistics() = default;

    private:
        std::mutex mutex_;
        std::vector<StateInfo> states_;
    };

// Retrieve current stats instance
    inline Statistics &stats() {
        return Statistics::instance();
    }

} // namespace floatTetWild
