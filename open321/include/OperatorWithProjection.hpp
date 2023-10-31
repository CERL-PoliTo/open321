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
* \file OperatorWithProjection.h
*/

#ifndef INCLUDE_OPERATORWITHPROJECTION_HPP
#define INCLUDE_OPERATORWITHPROJECTION_HPP

#include "open321_common.hpp"

namespace open321 {
   /**
    * @brief A template class for applying projection in combination with a discretized operator.
    *
    * The `OperatorWithProjection` class allows you to apply projection, defined by a `ProjectorT`,
    * in combination with a discretized operator, defined by an `OperatorT`. This class provides the
    * flexibility to choose between orthogonal or non-orthogonal projection based on the value of
    * the template boolean parameter `Orthogonal`.
    *
    * @tparam ProjectorT The type of the projector, which must implement the 'mvp' method.
    * @tparam OperatorT The type of the discretized operator, which must implement the 'mvp' method.
    * @tparam Orthogonal A boolean value (default: false) that determines whether orthogonal
    *                  projection is used. When set to true, the orthogonal projector is defined as
    *                   \f$P^* = I - P\f$, where \f$P$\f is the original projector.
    *
    * @note For the 'ProjectorT' and 'OperatorT' types, it is essential that they implement the 'mvp'
    * method, which typically represents a matrix-vector product.
    *
    * Example usage:
    * @code
    * // Create an instance of OperatorWithProjection using non-orthogonal projection.
    * OperatorWithProjection<NonOrthogonalProjector, DiscretizedOperator> op;
    *
    * // Apply the projection and operator to a vector 'x'.
    * Vector result = op.mvp(x);
    * @endcode
    */
    template<class ProjectorT, class OperatorT, bool Orthogonal = false>
    class OperatorWithProjection;
}

template<class ProjectorT, class OperatorT, bool Orthogonal>
struct Eigen::internal::traits<open321::OperatorWithProjection <ProjectorT, OperatorT, Orthogonal>> :
public Eigen::internal::traits<Eigen::SparseMatrix<open321::complex_t>> {};

namespace open321 {
    template<class ProjectorT, class OperatorT, bool Orthogonal>
    class OperatorWithProjection : public Eigen::EigenBase<OperatorWithProjection<ProjectorT, OperatorT, Orthogonal>> {
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

        ProjectorT &m_projector_matrix;
        OperatorT &m_operator_matrix;
        Eigen::Index m_n_bf;
        mutable VTypeC m_dest_Buf;

        Eigen::Index rows() const { return m_n_bf; }

        Eigen::Index cols() const { return m_n_bf; }

        VTypeC& getDestBuffer() const { return m_dest_Buf; }

        template<typename Rhs>
        Eigen::Product<OperatorWithProjection, Rhs, Eigen::AliasFreeProduct>
        operator*(const Eigen::MatrixBase<Rhs> &x) const {
            return Eigen::Product<OperatorWithProjection, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
        }

        OperatorWithProjection(ProjectorT &projector_matrix, OperatorT &operator_matrix) :
                m_projector_matrix(projector_matrix), m_operator_matrix(operator_matrix),
                m_n_bf(m_operator_matrix.cols()), m_dest_Buf(m_operator_matrix.cols()) {
        }

        template<typename Src, typename Dest, bool O = Orthogonal, typename std::enable_if<O>::type * = nullptr>
        void mvp(const Src &x, Dest &res) const {
            VTypeC resProjection = VTypeC::Zero(m_n_bf);
            m_projector_matrix.mvp(x, resProjection);
            VTypeC resOrthogonalProjection = x - resProjection;
            VTypeC resOpetatorProjectionRight = m_operator_matrix * resOrthogonalProjection;
            VTypeC res_tmp = VTypeC::Zero(m_n_bf);
            m_projector_matrix.mvp(resOpetatorProjectionRight, res_tmp);
            res = resOpetatorProjectionRight - res_tmp;
        }

        template<typename Src, typename Dest, bool O = Orthogonal, typename std::enable_if<!O>::type * = nullptr>
        void mvp(const Src &x, Dest &res) const {
            VTypeC resPSigma = VTypeC::Zero(m_n_bf);
            m_projector_matrix.mvp(x, resPSigma);
            VTypeC resTsPLambdaH = m_operator_matrix * resPSigma;
            VTypeC res_tmp = VTypeC::Zero(m_n_bf);
            m_projector_matrix.mvp(res, res_tmp);
        }
    };
}

namespace Eigen {
    namespace internal {
        template<class A, class B, bool C, typename Rhs>
        struct generic_product_impl<open321::OperatorWithProjection <
                                    A, B, C>, Rhs, SparseShape, DenseShape, GemvProduct>
        : generic_product_impl_base<open321::OperatorWithProjection < A, B, C>, Rhs, generic_product_impl<open321::OperatorWithProjection < A, B, C>, Rhs> > {
        typedef typename Product<open321::OperatorWithProjection < A, B, C>, Rhs>::Scalar Scalar;
            template<typename Dest>
            static void scaleAndAddTo(Dest &dst, const open321::OperatorWithProjection <A, B, C> &lhs, const Rhs &rhs, const Scalar &alpha) {
                open321::VTypeC dst_buffer = lhs.getDestBuffer();
                lhs.mvp(rhs, dst_buffer);
                dst.noalias() += alpha*dst_buffer;
            }
        };
    }
}

#endif
