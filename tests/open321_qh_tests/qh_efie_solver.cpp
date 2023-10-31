#include <iostream>

#include "constants.hpp"
#include "open321_common.hpp"
#include "QuasiHelmholtzEFIESolver.hpp"

#include <unsupported/Eigen/SparseExtra>

#include "qh_test_config.hpp"

using namespace open321;

TEST(open321Preconditioner, QuasiHelmholtzEFIESolver){
    floating_t kwave = 2 * pi * 3e8 / c0;

    std::string test_root(OPEN321_TEST_DATA_ROOT);

    std::string filename_th(test_root + "/Th.mtx");
    std::string filename_ts(test_root + "/Ts.mtx");
    std::string filename_rhs(test_root + "/rhs.mtx");
    std::string filename_rhs_extracted(test_root + "/rhs_extracted.mtx");
    std::string filename_sigma(test_root + "/Sigma.mtx");
    std::string filename_refsolution(test_root + "/solution.mtx");

    MtypeC th;
    MtypeC ts;
    MtypeC debug_m;
    VTypeC rhs;
    VTypeC rhs_extracted;
    MtypeSparse sigma;
    VTypeC ref_solution;
    Eigen::loadMarketDense(th, filename_th);
    Eigen::loadMarketDense(ts, filename_ts);
    Eigen::loadMarketVector(rhs, filename_rhs);
    Eigen::loadMarketVector(rhs_extracted, filename_rhs_extracted);
    Eigen::loadMarket(sigma, filename_sigma);
    Eigen::loadMarketVector(ref_solution, filename_refsolution);

    floating_t tolerance = 1e-5;
    QuasiHelmholtzEFIESolver qh_efie_solver(sigma, th, ts, kwave, tolerance, 200);
    VTypeC x = qh_efie_solver.solve<VTypeC>(rhs, rhs_extracted);

    EXPECT_TRUE(((x - ref_solution).norm() / ref_solution.norm() < tolerance * 3));
}
