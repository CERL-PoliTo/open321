// Copyright (C) 2023 Francesco P. Andriulli, Johann P. Bourhis, Damiano Franzò, Adrien Merlini
// All rights reserved.
//
// This file is part of Open321.
//
// Open321 is free software: you can redistribute it and/or modify it
// under the terms of the GNU Affero General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// Open321 is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
// License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with Open321.  If not, see <https://www.gnu.org/licenses/>.

/**
 * \file open321_common.h
 *
 * @brief Common Aliases for Eigen Library
 *
 * This header file defines various aliases and type definitions commonly used
 * throughout the project with the Eigen library. These aliases primarily refer to
 * Eigen's matrix and vector types. Using these aliases improves code readability
 * and maintainability, as they provide meaningful names for Eigen data types and
 * template parameters.
 *
 * This file is part of the project's core library and should be included wherever
 * Eigen matrix and vector types are used.
 */

#ifndef INCLUDE_OPEN321_COMMON_HPP
#define INCLUDE_OPEN321_COMMON_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "open321_config.hpp"

#include <cstddef>
#include <vector>
#include <iostream>
#include <complex>

namespace open321 {

    using floating_t = double;
    using complex_t = std::complex<floating_t>;

    constexpr complex_t J{0, 1};

    using Mtype = Eigen::Matrix<floating_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MtypeCol = Eigen::Matrix<floating_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using MtypeSparse = Eigen::SparseMatrix<floating_t, Eigen::RowMajor>;
    using MtypeSparseCol = Eigen::SparseMatrix<floating_t, Eigen::ColMajor>;
    using VType = Eigen::Matrix<floating_t, 1,Eigen::Dynamic, Eigen::RowMajor>;
    using VTypeCol = Eigen::Matrix<floating_t, Eigen::Dynamic, 1, Eigen::ColMajor>;

    using MtypeC = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MtypeCCol = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using MtypeSparseC = Eigen::SparseMatrix<complex_t, Eigen::RowMajor>;
    using MtypeSparseCCol = Eigen::SparseMatrix<complex_t, Eigen::ColMajor>;
    using VTypeC = Eigen::Matrix<complex_t, Eigen::Dynamic, 1, Eigen::ColMajor>;
    using VTypeCRow = Eigen::Matrix<complex_t, 1, Eigen::Dynamic, Eigen::RowMajor>;

}
#endif
