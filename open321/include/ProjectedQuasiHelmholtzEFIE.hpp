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
* \file ProjectedQuasiHelmholtzEFIE.h
*/

#ifndef INCLUDE_PROJECTEDQUASIHELMHOLTZEFIE_HPP
#define INCLUDE_PROJECTEDQUASIHELMHOLTZEFIE_HPP

#include "open321_common.hpp"
#include "QuasiHelmholtzProjector.hpp"
#include "OperatorWithProjection.hpp"
#include "utility.hpp"

namespace open321 {
    /**
     * \brief ProjectedQuasiHelmholtzEFIE - Projected Electric Field Integral Equation (EFIE)
     *
     * The `ProjectedQuasiHelmholtzEFIE` class is responsible for implementing the Electric Field Integral Equation (EFIE),
     * where quasi-Helmholtz projectors are applied.
     * @tparam THType  The type representing the matrix Th, which must implement the 'mvp' method for matrix-vector products.
     * @tparam TSType  The type representing the matrix Ts, which must implement the 'mvp' method for matrix-vector products.
     */
    template<class THType, class TSType>
    class ProjectedQuasiHelmholtzEFIE;
}

template<class THType, class TSType>
struct Eigen::internal::traits<open321::ProjectedQuasiHelmholtzEFIE<THType, TSType>> :
public Eigen::internal::traits<Eigen::SparseMatrix<open321::complex_t>> {};

namespace open321 {
    template<class THType, class TSType>
    class ProjectedQuasiHelmholtzEFIE : public Eigen::EigenBase<ProjectedQuasiHelmholtzEFIE<THType, TSType>> {
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

        THType &m_th;
        TSType &m_ts;
        QuasiHelmholtzProjector &m_psigma;
        size_t m_n_bf;
        floating_t m_kwave, m_c;
        mutable VTypeC m_dest_Buf;


        Eigen::Index rows() const { return m_n_bf; }

        Eigen::Index cols() const { return m_n_bf; }

        template<typename Rhs>
        Eigen::Product<ProjectedQuasiHelmholtzEFIE, Rhs, Eigen::AliasFreeProduct>
        operator*(const Eigen::MatrixBase<Rhs> &x) const {
            return Eigen::Product<ProjectedQuasiHelmholtzEFIE, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
        }

        ProjectedQuasiHelmholtzEFIE(THType &th, TSType &ts, QuasiHelmholtzProjector &psigma, floating_t kwave,
                                    floating_t c) :
                m_th(th), m_ts(ts), m_psigma(psigma), m_n_bf(th.cols()), m_kwave(kwave), m_c(c), m_dest_Buf(th.cols()) {
        }

        ProjectedQuasiHelmholtzEFIE(THType &th, TSType &ts, QuasiHelmholtzProjector &psigma, floating_t kwave) :
                m_th(th), m_ts(ts), m_psigma(psigma), m_n_bf(th.cols()), m_kwave(kwave) {
            OperatorWithProjection<QuasiHelmholtzProjector, MtypeC, true> pl_ts_pl(psigma, ts);
            floating_t num = power_method_highest_sv_symmetric_complex(th, 1e-2);
            floating_t denum = power_method_highest_sv_symmetric_complex(pl_ts_pl, 1e-2);
            m_c = std::sqrt(num / denum);
        }

        floating_t getC() {
            return m_c;
        }


        VTypeC& getDestBuffer() const
        {
            return m_dest_Buf;
        }

        // P_k = √(C/k) P_ΛH + i√(k/C) P_∑
        // T_k = ik T_s + 1/(ik) T_h

    /**
     * @brief Matrix-Vector Product (MVP) Operation
     *
     * This function performs a matrix-vector product operation where various projectors and scaling factors
     * are applied to the input vector 'x' to generate multiple result vectors 'res1' through 'res5'. It produces
     * the result of the operation the operation: \f$P_k T_k P_k x\f\f$, where:
     *
     * - \f$P_k = \sqrt{\frac{C}{k}} P_{\Lambda H} + i\sqrt{\frac{k}{C}} P_{\Sigma}\f$
     * - \f$T_k = i k T_s + \frac{1}{i k} T_h\f$
     *
     * The MVP operation is described as follows:
     *
     * - \f$res1 = P_{\Lambda H} T_s P_{\Lambda H} x\f, with scaling factor \f$C i\f$
     * - \f$res2 = P_{\Sigma} T_s P_{\Lambda H} x\f, with scaling factor \f$-k\f$
     * - \f$res3 = P_{\Lambda H} T_s P_{\Sigma} x\f, with scaling factor \f$-k\f$
     * - \f$res4 = P_{\Sigma} T_s P_{\Sigma} x\f, with scaling factor \f$-i k^2 / C\f$
     * - \f$res5 = T_h x\f, with scaling factor \f$i / C\f$
     *
     * The result of the mvp is computed as the sum of all the vectors 'res1' through 'res5'
     *
     * @param x The input vector for the MVP operation.
     * @param res The destination vector for the results of the MVP operation.
     *
     * @tparam Src The source vector type.
     * @tparam Dest The destination vector type.
     */
    template<typename Src, typename Dest>
    void mvp(const Src &x, Dest &res) const {
        // P_∑ x
        // P_ΛH x = (1 - P_∑) x
        VTypeC resPSigma = VTypeC::Zero(m_n_bf);
        m_psigma.mvp(x, resPSigma);
        VTypeC resPLambdaH = x - resPSigma;

        //T_s P_ΛH
        //T_s P_∑
        VTypeC resTsPLambdaH = m_ts * resPLambdaH;
        VTypeC resTsPSigma = m_ts * resPSigma;

        VTypeC res2 = VTypeC::Zero(m_n_bf);
        m_psigma.mvp(resTsPLambdaH, res2);
        VTypeC res1 = resTsPLambdaH - res2;

        VTypeC res4 = VTypeC::Zero(m_n_bf);
        m_psigma.mvp(resTsPSigma, res4);
        VTypeC res3 = resTsPSigma - res4;

        VTypeC res5 = m_th * x;

        VTypeC resPartial = (m_c * J) * res1.array() - m_kwave * res2.array()
                            - m_kwave * res3.array() - (J * m_kwave * m_kwave / m_c) * res4.array()
                            + (J / m_c) * res5.array();
        res = resPartial;
    }

    };
}

namespace Eigen {
    namespace internal {
        template<class A, class B, typename Rhs>
        struct generic_product_impl<open321::ProjectedQuasiHelmholtzEFIE<A, B>, Rhs, SparseShape, DenseShape, GemvProduct>
                : generic_product_impl_base<open321::ProjectedQuasiHelmholtzEFIE<A, B>, Rhs, generic_product_impl<open321::ProjectedQuasiHelmholtzEFIE< A, B>, Rhs> > {
            typedef typename Product<open321::ProjectedQuasiHelmholtzEFIE<A, B>, Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest &dst, const open321::ProjectedQuasiHelmholtzEFIE<A, B> &lhs, const Rhs &rhs,
                                      const Scalar &alpha) {
                open321::VTypeC dst_buffer = lhs.getDestBuffer();
                lhs.mvp(rhs, dst_buffer);
                dst.noalias() += alpha*dst_buffer;
            }
        };
    }
}

#endif
