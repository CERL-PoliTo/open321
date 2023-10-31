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
* \file utility.h
* Contains some utility functions
*/
#ifndef INCLUDE_UTILITY_HPP
#define INCLUDE_UTILITY_HPP

#include "open321_common.hpp"

namespace open321 {
    /**
     *
     * @brief Power Method for Computing the Highest Singular Value of a Symmetric Complex Matrix
     *
     * The `power_method_highest_sv_symmetric_complex` method employs the power iteration method to compute
     * the highest singular value of a given symmetric complex matrix 'm'.
     *
     * @param m The symmetric complex matrix for which the highest singular value is to be computed.
     * @param tolerance The convergence tolerance for stopping the power iteration.
     * @param num_iterations (Optional) The maximum number of iterations (default: 100).
     *
     * @return The highest singular value of the matrix 'm'.
     *
     * This method iteratively applies the power iteration technique to estimate the highest singular value of the
     * given symmetric complex matrix 'm'. It terminates when the estimate converges within the specified 'tolerance'
     * or after the maximum number of iterations is reached.
     *
     * Example usage:
     * @code
     * // Create a symmetric complex matrix 'M' and specify the convergence tolerance.
     * SymmetricComplexMatrixType M;
     * floating_t tolerance = 1e-6;
     *
     * // Compute the highest singular value using power iteration.
     * floating_t highestSV = power_method_highest_sv_symmetric_complex(M, tolerance);
     * @endcode
     */
    template<typename M>
    floating_t power_method_highest_sv_symmetric_complex(const M &m, floating_t tolerance, size_t num_iterations=100) {
        VTypeC x = VTypeC::Random(m.cols());
        VTypeC x_next = VTypeC::Zero(m.cols());
        x.normalize();

        floating_t current_largest_sv, largest_sv = 1e100;
        bool convergence = false;

        for (size_t i = 0; (i < num_iterations) & (!convergence); i++) {
            // Apply the power iteration
            x_next.noalias() = m * x;
            current_largest_sv = x_next.norm();
            if (std::abs(current_largest_sv - largest_sv) < tolerance)
                convergence = true;
            largest_sv = current_largest_sv;
            // Normalize the vector at each iteration
            x_next.normalize();
            x = x_next;
        }

        // Calculate the largest singular value
        largest_sv = (m * x).norm();

        return largest_sv;
    }
}
#endif
