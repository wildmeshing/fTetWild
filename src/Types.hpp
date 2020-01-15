// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once


#include <Eigen/Dense>

#include <nlohmann/json.hpp>

namespace floatTetWild {
    //fortest
    extern double eps_zero;//1e-8 ~ 1e-16, *10
    extern double eps_prep_r;//0 ~ 0.9, +0.1
    extern int snap_shrink;//1e-0 ~ 1e-9, *10
    //fortest

#ifdef FLOAT_TETWILD_USE_FLOAT
    typedef float Scalar;
#define SCALAR_ZERO 1e-6
#define SCALAR_ZERO_2 1e-12
#define SCALAR_ZERO_3 1e-18
#else
    typedef double Scalar;
//#define SCALAR_ZERO 1e-8
//#define SCALAR_ZERO_2 1e-16
//#define SCALAR_ZERO_3 1e-24
#define SCALAR_ZERO eps_zero
#define SCALAR_ZERO_2 eps_zero*eps_zero
#define SCALAR_ZERO_3 eps_zero*eps_zero*eps_zero
#endif

#define STORE_SAMPLE_POINTS

    // Json
    using json = nlohmann::json;

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;

    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Matrix<Scalar, 2, 1> Vector2;


    typedef Eigen::Matrix<int, 4, 1> Vector4i;
    typedef Eigen::Matrix<int, 3, 1> Vector3i;
    typedef Eigen::Matrix<int, 2, 1> Vector2i;
}
