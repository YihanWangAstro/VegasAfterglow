#include <boost/test/unit_test.hpp>
#include <cmath>

#include "dynamics/reverse-shock.hpp"
#include "util/macros.h"

// Integration-level tests of the coupled forward-reverse shock system
// (FRShockEqn): solve the ODEs through crossing and deceleration and assert
// physical invariants no pointwise unit test can see. These are the tests
// that would have caught the 2026-07 magnetized-shell runaway (NaN at the
// weak-penetration singularity silently freezing dGamma/dt).
//
// Energy-budget calibration (adiabatic limit, eps_e = 1e-6): the blast wave
// ends with E_blast ~= 1.13 * eps_k / (1 + sigma) -- the kinetic share of the
// shell energy plus swept-mass bookkeeping; the magnetic fraction
// sigma/(1+sigma) is currently NOT transferred to the blast (documented
// behavior, see the method paper Eqs. 39-44 discussion).

BOOST_AUTO_TEST_SUITE(ReverseShockIntegration)

namespace {

    struct PairRunResult {
        bool all_finite = true;
        bool crossed = false;
        int gamma_rise_violations = 0; // post-crossing Gamma increases
        Real Gamma_final = 0;
        Real m3_over_m4 = 0;
        Real E_blast_over_eps_k = 0;
    };

    PairRunResult run_pair_system(Real sigma0) {
        const Real eps_k = Real(1e53) * unit::erg / (4 * con::pi);
        Ejecta jet([eps_k](Real, Real) { return eps_k; }, [](Real, Real) { return Real(300.0); },
                   [sigma0](Real, Real) { return sigma0; }, func::zero_3d, func::zero_3d, false, 1 * unit::sec);
        ISM medium(1.0 / unit::cm3);
        RadParams rad;
        rad.eps_e = 1e-6; // adiabatic limit: negligible radiative losses
        rad.eps_B = 0.01;
        rad.p = 2.3;
        rad.xi_e = 1.0;

        FRShockEqn eqn(medium, jet, 0.0, 0.01, rad, rad);
        typename decltype(eqn)::State state;
        const Real t0 = 1e-5 * unit::sec;
        eqn.set_init_state(state, t0);

        using namespace boost::numeric::odeint;
        auto stepper = make_dense_output(1e-6, 1e-6, runge_kutta_dopri5<typename decltype(eqn)::State>());
        stepper.initialize(state, t0, 1e-9 * t0);

        PairRunResult res;
        const Real t_end = 1e8 * unit::sec;
        Real Gamma_prev = 0;

        for (size_t steps = 0; stepper.current_time() <= t_end && steps < 200000; ++steps) {
            stepper.do_step(eqn);
            auto const& s = stepper.current_state();
            if (!std::isfinite(s.Gamma + s.m3 + s.U3_th + s.x3 + s.U2_th)) {
                res.all_finite = false;
                break;
            }
            if (!res.crossed && s.m3 >= 0.999 * s.m4 && stepper.current_time() > 1.5 * unit::sec) {
                res.crossed = true;
                Gamma_prev = s.Gamma;
            }
            if (res.crossed) {
                if (s.Gamma > Gamma_prev * (1 + 1e-10)) {
                    ++res.gamma_rise_violations;
                }
                Gamma_prev = s.Gamma;
            }
        }

        auto const& s = stepper.current_state();
        res.Gamma_final = s.Gamma;
        res.m3_over_m4 = (s.m4 > 0) ? s.m3 / s.m4 : 0;
        const Real ad2 = physics::thermo::adiabatic_idx(s.Gamma);
        const Real Geff = compute_effective_Gamma(ad2, s.Gamma);
        const Real E_blast = (s.Gamma - 1) * (s.m2 + s.m3) * con::c2 + Geff * (s.U2_th + s.U3_th);
        res.E_blast_over_eps_k = E_blast / eps_k;
        return res;
    }

} // namespace

// ---------------------------------------------------------------------------
// Coupled system stays finite, crosses, decelerates monotonically, and lands
// on the expected energy budget -- for unmagnetized and magnetized shells.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(pair_system_physical_invariants) {
    for (Real sigma : {0.0, 0.1, 1.0}) {
        auto res = run_pair_system(sigma);
        BOOST_TEST_CONTEXT("sigma = " << sigma) {
            BOOST_CHECK(res.all_finite);
            BOOST_CHECK(res.crossed);
            BOOST_CHECK_CLOSE(res.m3_over_m4, 1.0, 0.2); // crossing completes
            BOOST_CHECK_EQUAL(res.gamma_rise_violations, 0);
            BOOST_CHECK_LT(res.Gamma_final, 1.2); // fully decelerated by t_end
            // Energy budget: kinetic share of the shell energy reaches the blast
            // (calibrated 1.13/(1+sigma); band allows integration/EOS spread)
            const Real scaled = res.E_blast_over_eps_k * (1 + sigma);
            BOOST_CHECK_GT(scaled, 1.0);
            BOOST_CHECK_LT(scaled, 1.3);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
