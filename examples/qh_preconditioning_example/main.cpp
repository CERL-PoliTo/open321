#include <iostream>

#include "open321_common.hpp"
#include "constants.hpp"
#include "QuasiHelmholtzProjector.hpp"
#include "ProjectedQuasiHelmholtzEFIE.hpp"
#include "OperatorWithProjection.hpp"
#include "utility.hpp"

#include <unsupported/Eigen/SparseExtra>

#include "qh_example_config.hpp"

using namespace open321;

int main() {
    floating_t kwave = 2 * pi * 3e8 / c0;

    std::string example_root(OPEN321_EXAMPLE_DATA_ROOT);

    std::string filename_th(example_root+"/Th.mtx");
    std::string filename_ts(example_root+"/Ts.mtx");
    std::string filename_rhs(example_root+"/rhs.mtx");
    std::string filename_rhs_extracted(example_root+"/rhs_extracted.mtx");
    std::string filename_sigma(example_root+"/Sigma.mtx");
    std::string filename_refsolution(example_root+"/solution.mtx");

    MtypeC th;
    MtypeC ts;
    MtypeC debug_m;
    VTypeC rhs;
    VTypeC rhs_extracted;
    MtypeSparse sigma;
    VTypeC solution;
    Eigen::loadMarketDense(th, filename_th);
    Eigen::loadMarketDense(ts, filename_ts);
    Eigen::loadMarketVector(rhs, filename_rhs);
    Eigen::loadMarketVector(rhs_extracted, filename_rhs_extracted);
    Eigen::loadMarket(sigma, filename_sigma);
    Eigen::loadMarketVector(solution, filename_refsolution);

    Eigen::Index nbf = th.rows();

    // P_k = √(C/k) P_ΛH + i √(k/C) P_∑
    // T_k = ik T_s + 1/(ik) T_h
    //Estimate Build projector class

    floating_t tolerance = 1e-5;
    QuasiHelmholtzProjector psigma(sigma, tolerance);

    //Estimate C
    OperatorWithProjection<QuasiHelmholtzProjector, MtypeC, true> pl_ts_pl(psigma, ts);
    floating_t num = power_method_highest_sv_symmetric_complex(th, tolerance);
    floating_t denum = power_method_highest_sv_symmetric_complex(pl_ts_pl, tolerance);

    std:: cout << "C num : " << num <<std::endl;
    std:: cout << "C denum : " << denum <<std::endl;

    floating_t c = sqrt(num / denum);

    // RHS
    VTypeC rhs_psigma = VTypeC::Zero(nbf);
    psigma.mvp(rhs, rhs_psigma);

    VTypeC rhs_ext_psigma = VTypeC::Zero(nbf);
    psigma.mvp(rhs_extracted, rhs_ext_psigma);

    VTypeC rhs_all = sqrt(c/kwave) * (rhs_extracted - rhs_ext_psigma).array() + rhs_psigma.array() * (sqrt(kwave / c)  * J);

    ProjectedQuasiHelmholtzEFIE qh(th, ts, psigma, kwave, c);

    Eigen::GMRES<ProjectedQuasiHelmholtzEFIE<MtypeC, MtypeC>, Eigen::IdentityPreconditioner> gmres_solver(qh);
    gmres_solver.set_restart(200);
    gmres_solver.setTolerance(tolerance);
    VTypeC x = VTypeC::Zero(nbf);
    x = gmres_solver.solve(rhs_all);

    VTypeC psigma_x = VTypeC::Zero(nbf);
    psigma.mvp(x, psigma_x);
    VTypeC y_star = J*sqrt(kwave/c) * psigma_x;
    VTypeC y_loop = sqrt(c/kwave) * (x - psigma_x);
    VTypeC y = y_star+y_loop;

    floating_t rel_err = (y-solution).norm() / solution.norm();
    floating_t abs_err = (y-solution).norm();
    std::cout << "GMRES: Iterations : " << gmres_solver.iterations() << std::endl;
    std::cout << "GMRES: rerror : " << gmres_solver.error() << std::endl;
    std::cout << rel_err << std::endl;
    std::cout << abs_err << std::endl;

    return 0;
}
