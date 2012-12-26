#include <src/rysint/_gvrr_drv.h>
void bagel::GVRR_Driver::gvrr_driver_2_1_6_4(double* out, const double* const roots, const double* const weights, const double& coeff,
    const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
    const double* const p, const double* const q, const double& xp, const double& xq, const size_t& size_block,
    const double* const expo, const double* const transx, const double* const transy, const double* const transz,
    const double* const trans2x, const double* const trans2y, const double* const trans2z, double* const intermediate,
    double* const final_x, double* const final_y, double* const final_z,
    double* const final_xa, double* const final_xb, double* const final_xc,
    double* const final_ya, double* const final_yb, double* const final_yc,
    double* const final_za, double* const final_zb, double* const final_zc,
    double* const workx, double* const worky, double* const workz) {
  bagel::gvrr_driver<2,1,6,4,8>(out, roots, weights, coeff, a, b, c, d, p, q, xp, xq, size_block,
    expo, transx, transy, transz, trans2x, trans2y, trans2z, intermediate, final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc,
    final_za, final_zb, final_zc, workx, worky, workz);
}
