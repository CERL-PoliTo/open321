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
* \file QuasiHelmholtzProjector.h
*/
#ifndef INCLUDE_QUASIHELMHOLTZPROJECTOR_HPP
#define INCLUDE_QUASIHELMHOLTZPROJECTOR_HPP

#include <unsupported/Eigen/IterativeSolvers>

#include "open321_common.hpp"

#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>

#include <amgcl/backend/eigen.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/value_type/eigen.hpp>
#include <amgcl/value_type/complex.hpp>

namespace open321 {
    typedef amgcl::backend::eigen<complex_t> Backend;
    typedef amgcl::make_solver<
    amgcl::amg<
            Backend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
    >,
    amgcl::solver::cg<Backend>
    > Solver;

    /**
     * @brief Quasi-Helmholtz Projector
     *
     * The `QuasiHelmholtzProjector` class implements a Quasi-Helmholtz projector for a given matrix,
     * such as 'Sigma'. This projector is designed to improve numerical stability in electromagnetic simulations
     * and can be used in combination with iterative solvers.
     *
     * The class inherits traits from `Eigen::SparseMatrix` to ensure compatibility with iterative solvers.
     *
     * @tparam MatrixType Type of the matrix for which the QH projector is implemented.
     *
     * Example usage:
     * @code
     * // Create an instance of QuasiHelmholtzProjector for a specific matrix, e.g., 'Sigma'.
     * QuasiHelmholtzProjector<SigmaMatrixType> qhProjector;
     *
     * // Apply the QH projector to a vector 'x'.
     * Eigen::VectorXd projectedVector = qhProjector * x;
     * @endcode
     */
    class QuasiHelmholtzProjector;
}

template<>
struct Eigen::internal::traits<open321::QuasiHelmholtzProjector> :
        public Eigen::internal::traits<Eigen::SparseMatrix<open321::floating_t>> {};

namespace open321 {
    class QuasiHelmholtzProjector : public Eigen::EigenBase<QuasiHelmholtzProjector> {

        Solver::params settings_prm(floating_t tolerance) {
            Solver::params prm;
            prm.solver.tol = tolerance;
            prm.solver.maxiter = 10;
            return prm;
        }

    public:
        // Required typedefs, constants, and method:
        typedef complex_t Scalar;
        typedef double RealScalar;
        typedef int StorageIndex;
        enum {
            ColsAtCompileTime = Eigen::Dynamic,
            MaxColsAtCompileTime = Eigen::Dynamic,
            IsRowMajor = false
        };

        MtypeSparse m_Sigma;
        Eigen::Index n_bf;
        floating_t m_tolerance;
        MtypeSparseC m_graphlaplacian; //tmp
        Solver solve;

        mutable VTypeC prod, rhs, dst_buffer;

        Eigen::Index rows() const { return n_bf; }

        Eigen::Index cols() const { return n_bf; }

        VTypeC& getBuffer() const {
            return dst_buffer;
        }

        template<typename Rhs>
        Eigen::Product<QuasiHelmholtzProjector, Rhs, Eigen::AliasFreeProduct>
        operator*(const Eigen::MatrixBase<Rhs> &x) const {
            return Eigen::Product<QuasiHelmholtzProjector, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
        }

        QuasiHelmholtzProjector(open321::MtypeSparse &Sigma, open321::floating_t tolerance) :
                m_Sigma(Sigma), n_bf(Sigma.rows()), m_tolerance(tolerance),
                m_graphlaplacian((Sigma.transpose() * Sigma).pruned().cast<complex_t>()),
                solve(m_graphlaplacian, settings_prm(tolerance)), prod(Sigma.cols()),rhs(Sigma.cols()) {
        }

        template<typename X, typename Dest>
        void mvp(const X &x, Dest &res) const {
            prod.setZero();
            rhs.noalias() = m_Sigma.transpose() * x; //TODO same
            solve(rhs, prod);
            res.noalias() = m_Sigma * prod;
        }

    };
}

namespace Eigen {
    namespace internal {
        template<typename Rhs>
        struct generic_product_impl<open321::QuasiHelmholtzProjector, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
                : generic_product_impl_base<open321::QuasiHelmholtzProjector, Rhs, generic_product_impl<open321::QuasiHelmholtzProjector, Rhs> > {
            typedef typename Product<open321::QuasiHelmholtzProjector, Rhs>::Scalar Scalar;

            template<typename Dest>
            static void
            scaleAndAddTo(Dest &dst, const open321::QuasiHelmholtzProjector &lhs, const Rhs &rhs, const Scalar &alpha) {
                open321::VTypeC& dst_buffer = lhs.getBuffer();
                lhs.mvp(rhs, dst_buffer);
                dst.noalias() += alpha * dst_buffer;
            }
        };
    }
}
#endif
