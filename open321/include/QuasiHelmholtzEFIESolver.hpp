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
* \file QuasiHelmholtzEFIESolver.h
*/
#ifndef INCLUDE_QUASIHELMHOLTZEFIESOLVER_HPP
#define INCLUDE_QUASIHELMHOLTZEFIESOLVER_HPP

#include <unsupported/Eigen/IterativeSolvers>

#include "open321_common.hpp"
#include "ProjectedQuasiHelmholtzEFIE.hpp"
#include "QuasiHelmholtzProjector.hpp"


namespace open321 {
    /**
     * @brief Solver for Quasi-Helmholtz Electric Field Integral Equation (EFIE)
     *
     *
     * The `QuasiHelmholtzEFIESolver` class provides a solver for the Quasi-Helmholtz Electric Field Integral Equation (EFIE).
     * It utilizes `Eigen::GMRES` to solve the problem, and its primary method is 'solve'.
     *
     * @tparam THType Type representing the matrix Th, which must implement the 'mvp' method.
     * @tparam TSType Type representing the matrix Ts, which must implement the 'mvp' method.
     *
     * The `solve` method can be used in two ways:
     * 1. `solve(rhs)` to solve the EFIE equation using the right-hand side vector 'rhs'.
     * 2. `solve(rhs, rhs_extracted)` to solve the EFIE equation using 'rhs' and store the extracted solution in 'rhs_extracted'.
     *
     * Example usage:
     * @code
     * // Create a QuasiHelmholtzEFIESolver instance with default Eigen::GMRES solver.
     * QuasiHelmholtzEFIESolver<ThMatrixType, TsMatrixType> solver;
     * // Define right-hand side vector and solution vector.
     * Eigen::VectorXd rhs, solution;
     * // Solve the EFIE using 'solve' method.
     * solver.solve(rhs, solution);
     * @endcode
     */
    template<class THType, class TSType>
    class QuasiHelmholtzEFIESolver;

    template<class THType, class TSType>
    class QuasiHelmholtzEFIESolver {
    public:
        QuasiHelmholtzProjector psigma;
        ProjectedQuasiHelmholtzEFIE<THType, TSType> m_precond_efie;
        Eigen::Index m_n;
        floating_t m_c, m_tolerance, m_kwave;
        size_t m_n_restart;

        QuasiHelmholtzEFIESolver(MtypeSparse &sigma_matrix, THType &th, TSType &ts,
                                 floating_t kwave, floating_t tolerance, size_t n_restart=30) :
                psigma(sigma_matrix, tolerance),
                m_precond_efie(th, ts, psigma, kwave), m_n(m_precond_efie.cols()),
                m_c(m_precond_efie.getC()), m_tolerance(tolerance), m_kwave(kwave), m_n_restart(n_restart) {
        }

        template<typename DestType, typename RHSType>
        DestType solve(const RHSType &rhs) {
            RHSType rhs_psigma = RHSType::Zero(m_n);
            psigma.mvp(rhs, rhs_psigma);

            RHSType rhs_all = sqrt(m_c / m_kwave) * (rhs - rhs_psigma).array() + rhs_psigma.array() * (sqrt(m_c / m_kwave) * J);

            Eigen::GMRES<ProjectedQuasiHelmholtzEFIE<THType, TSType>, Eigen::IdentityPreconditioner> gmres_solver(m_precond_efie);
            gmres_solver.set_restart(m_n_restart);
            gmres_solver.setTolerance(m_tolerance);

            DestType x(m_n);
            x = gmres_solver.solve(rhs_all);

            DestType psigma_x = DestType::Zero(m_n);
            psigma.mvp(x, psigma_x);
            DestType y_star = J * sqrt(m_kwave / m_c) * psigma_x;
            DestType y_loop = sqrt(m_c / m_kwave) * (x - psigma_x);
            DestType y = y_star + y_loop;
            return y;

        }

        template<typename DestType, typename RHSType>
        DestType solve(const RHSType &rhs, const RHSType &rhs_extracted) {
            RHSType rhs_psigma = RHSType::Zero(m_n);
            psigma.mvp(rhs, rhs_psigma);

            RHSType rhs_ext_psigma = RHSType::Zero(m_n);
            psigma.mvp(rhs_extracted, rhs_ext_psigma);

            RHSType rhs_all = sqrt(m_c / m_kwave) * (rhs_extracted - rhs_ext_psigma).array() +
                              rhs_psigma.array() * (sqrt(m_kwave / m_c) * J);

            Eigen::GMRES<ProjectedQuasiHelmholtzEFIE<THType, TSType>, Eigen::IdentityPreconditioner> gmres_solver(m_precond_efie);
            gmres_solver.set_restart(m_n_restart);
            gmres_solver.setTolerance(m_tolerance);

            DestType x = DestType::Zero(m_n);
            x = gmres_solver.solve(rhs_all);

            DestType psigma_x = DestType::Zero(m_n);
            psigma.mvp(x, psigma_x);
            DestType y_star = J * sqrt(m_kwave / m_c) * psigma_x;
            DestType y_loop = sqrt(m_c / m_kwave) * (x - psigma_x);
            DestType y = y_star + y_loop;
            return y;
        }
    };

}
#endif
