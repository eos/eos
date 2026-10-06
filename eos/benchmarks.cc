/*
 * Copyright (c) 2026 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/benchmarks.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/observable-impl.hh>
#include <eos/signal-pdf-impl.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/parameters.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace eos
{
    namespace
    {
        /*
         * Benchmark B1 in x = m_mumu: a constant amplitude interfering with two narrow poles within a box,
         *
         *   P(x) = |1 + R_1 / (x - z_1) + R_2 / (x - z_2)|^2  for 0 <= x <= W, and zero otherwise,
         *
         * with z_j = x_j - i gamma_j. Near z_j, P has the singular term Re[c_j / (x - z_j)], with
         * c_j = 2 R_j conj(A_j) + i |R_j|^2 / gamma_j and A_j = 1 + R_k / (conj(z_j) - z_k), k != j.
         * The normalization is known in closed form.
         */
        class BenchmarkB1 : public ParameterUser
        {
            public:
                static const std::vector<OptionSpecification> options;

                static const std::set<ReferenceName> references;

                static const std::string description;

                UsedParameter x_1, gamma_1, x_2, gamma_2;

                UsedParameter abs_R_1, arg_R_1, abs_R_2, arg_R_2;

                UsedParameter W;

                BenchmarkB1(const Parameters & p, const Options &) :
                    x_1(p["B1::x_1@RvDV:2026A"], *this),
                    gamma_1(p["B1::gamma_1@RvDV:2026A"], *this),
                    x_2(p["B1::x_2@RvDV:2026A"], *this),
                    gamma_2(p["B1::gamma_2@RvDV:2026A"], *this),
                    abs_R_1(p["B1::abs{R_1}@RvDV:2026A"], *this),
                    arg_R_1(p["B1::arg{R_1}@RvDV:2026A"], *this),
                    abs_R_2(p["B1::abs{R_2}@RvDV:2026A"], *this),
                    arg_R_2(p["B1::arg{R_2}@RvDV:2026A"], *this),
                    W(p["B1::W@RvDV:2026A"], *this)
                {
                }

                std::array<complex<double>, 2>
                residues() const
                {
                    return { std::polar(abs_R_1(), arg_R_1()), std::polar(abs_R_2(), arg_R_2()) };
                }

                std::array<complex<double>, 2>
                poles() const
                {
                    return { complex<double>(x_1(), -gamma_1()), complex<double>(x_2(), -gamma_2()) };
                }

                // the coefficient c_j of the singular term of P at z_j
                complex<double>
                coefficient(const std::size_t & j) const
                {
                    const auto            R = residues();
                    const auto            z = poles();
                    const std::size_t     k = 1 - j;
                    const complex<double> A = 1.0 + R[k] / (std::conj(z[j]) - z[k]);

                    return 2.0 * R[j] * std::conj(A) + 1.0i * std::norm(R[j]) / (-z[j].imag());
                }

                double
                pdf(const double & x) const
                {
                    if ((x < 0.0) || (x > W()))
                    {
                        return 0.0;
                    }

                    const auto R = residues();
                    const auto z = poles();

                    return std::norm(1.0 + R[0] / (x - z[0]) + R[1] / (x - z[1]));
                }

                double
                normalization(const double & x_min, const double & x_max) const
                {
                    const double lo = std::max(x_min, 0.0), hi = std::min(x_max, W());
                    if (lo >= hi)
                    {
                        return 0.0;
                    }

                    const auto R = residues();
                    const auto z = poles();

                    // with p_j = R_j / (x - z_j): int p_j = R_j (L_j - i T_j), and int |p_j|^2 = |R_j|^2 T_j / gamma_j
                    std::array<double, 2> L, T;
                    for (std::size_t j = 0; j < 2; ++j)
                    {
                        const double x = z[j].real(), gamma = -z[j].imag();
                        L[j] = 0.5 * std::log((power_of<2>(hi - x) + gamma * gamma) / (power_of<2>(lo - x) + gamma * gamma));
                        T[j] = std::atan((hi - x) / gamma) - std::atan((lo - x) / gamma);
                    }

                    // int 2 Re[conj(p_1) p_2] = 2 Re[K ((L_2 - L_1) - i (T_1 + T_2))]
                    const complex<double> K = std::conj(R[0]) * R[1] / (z[1] - std::conj(z[0]));

                    double result = hi - lo;
                    for (std::size_t j = 0; j < 2; ++j)
                    {
                        result += 2.0 * (R[j] * complex<double>(L[j], -T[j])).real() + std::norm(R[j]) * T[j] / (-z[j].imag());
                    }

                    return result + 2.0 * (K * complex<double>(L[1] - L[0], -(T[0] + T[1]))).real();
                }

                double
                re_mu_1() const
                {
                    return x_1();
                }

                double
                im_mu_1() const
                {
                    return -gamma_1();
                }

                double
                re_mu_2() const
                {
                    return x_2();
                }

                double
                im_mu_2() const
                {
                    return -gamma_2();
                }

                double
                re_c_1() const
                {
                    return coefficient(0).real();
                }

                double
                im_c_1() const
                {
                    return coefficient(0).imag();
                }

                double
                re_c_2() const
                {
                    return coefficient(1).real();
                }

                double
                im_c_2() const
                {
                    return coefficient(1).imag();
                }

                static ObservableEntry::OptionIterator
                begin_options()
                {
                    return options.begin();
                }

                static ObservableEntry::OptionIterator
                end_options()
                {
                    return options.end();
                }
        };

        const std::vector<OptionSpecification> BenchmarkB1::options{};

        const std::set<ReferenceName> BenchmarkB1::references{};

        const std::string BenchmarkB1::description = "Benchmark B1: a constant amplitude interfering with two narrow poles within a box in m_mumu.";

        /*
         * Benchmark B2 in m = m_mumu and c = cos(theta_l), with q^2 = m^2:
         *
         *   P(m, c) = Phi(m) [ 3/4 |a_V|^2 (1 - c^2) + 1/2 |a_S + a_T c|^2 ],
         *   Phi(m)  = 2 m lambda^{3/2}(m_B^2, m_K^2, m^2) sqrt(1 - 4 m_mu^2 / m^2),
         *
         * within the phase space and |c| <= 1. As in B -> K l l, the scalar amplitude interferes with a tensor-like
         * amplitude through the lepton structure (a_S + a_T c), which keeps P non-negative for all parameters. The
         * resonances R_j / (q^2 - s_j) in a_V are poles in m at mu_j = sqrt(s_j), with residues rho_j = R_j / (2 mu_j).
         * Near mu_j, Phi |a_V|^2 has the singular term Re[c_j / (m - mu_j)], with
         * c_j = Phi(mu_j) [2 conj(a_V^{(j)}(conj(mu_j))) rho_j + i |rho_j|^2 / gamma_j],
         * gamma_j = -Im mu_j, and a_V^{(j)} = a_V - rho_j / (m - mu_j). The normalization integrates the
         * remainder numerically, and the singular terms in closed form.
         */
        class BenchmarkB2 : public ParameterUser
        {
            public:
                static const std::vector<OptionSpecification> options;

                static const std::set<ReferenceName> references;

                static const std::string description;

                UsedParameter m_B, m_K, m_mu, m_Bsst, m_Bs0;

                UsedParameter m_1, Gamma_1, m_2, Gamma_2;

                UsedParameter V_0, S_0, abs_T_0_over_S_0, arg_T_0;

                UsedParameter abs_R_1, arg_R_1, abs_R_2, arg_R_2;

                BenchmarkB2(const Parameters & p, const Options &) :
                    m_B(p["B2::m_B@RvDV:2026A"], *this),
                    m_K(p["B2::m_K@RvDV:2026A"], *this),
                    m_mu(p["B2::m_mu@RvDV:2026A"], *this),
                    m_Bsst(p["B2::m_B_s^*@RvDV:2026A"], *this),
                    m_Bs0(p["B2::m_B_s0@RvDV:2026A"], *this),
                    m_1(p["B2::m_1@RvDV:2026A"], *this),
                    Gamma_1(p["B2::Gamma_1@RvDV:2026A"], *this),
                    m_2(p["B2::m_2@RvDV:2026A"], *this),
                    Gamma_2(p["B2::Gamma_2@RvDV:2026A"], *this),
                    V_0(p["B2::V_0@RvDV:2026A"], *this),
                    S_0(p["B2::S_0@RvDV:2026A"], *this),
                    abs_T_0_over_S_0(p["B2::abs{T_0}/S_0@RvDV:2026A"], *this),
                    arg_T_0(p["B2::arg{T_0}@RvDV:2026A"], *this),
                    abs_R_1(p["B2::abs{R_1}@RvDV:2026A"], *this),
                    arg_R_1(p["B2::arg{R_1}@RvDV:2026A"], *this),
                    abs_R_2(p["B2::abs{R_2}@RvDV:2026A"], *this),
                    arg_R_2(p["B2::arg{R_2}@RvDV:2026A"], *this)
                {
                }

                double
                m_min() const
                {
                    return 2.0 * m_mu();
                }

                double
                m_max() const
                {
                    return m_B() - m_K();
                }

                complex<double>
                phase_space(const complex<double> & m) const
                {
                    const complex<double> q2 = m * m;
                    const complex<double> lambda =
                            power_of<2>(m_B() * m_B()) + power_of<2>(m_K() * m_K()) + q2 * q2 - 2.0 * (m_B() * m_B() * m_K() * m_K() + m_B() * m_B() * q2 + m_K() * m_K() * q2);

                    return 2.0 * m * std::pow(lambda, 1.5) * std::sqrt(1.0 - 4.0 * m_mu() * m_mu() / q2);
                }

                std::array<complex<double>, 2>
                residues() const
                {
                    return { std::polar(abs_R_1(), arg_R_1()), std::polar(abs_R_2(), arg_R_2()) };
                }

                std::array<complex<double>, 2>
                poles() const
                {
                    return { std::sqrt(complex<double>(m_1() * m_1(), -m_1() * Gamma_1())), std::sqrt(complex<double>(m_2() * m_2(), -m_2() * Gamma_2())) };
                }

                complex<double>
                a_V(const complex<double> & m) const
                {
                    const complex<double> q2 = m * m;
                    const auto            R  = residues();

                    return V_0() / (1.0 - q2 / (m_Bsst() * m_Bsst())) + R[0] / (q2 - complex<double>(m_1() * m_1(), -m_1() * Gamma_1()))
                           + R[1] / (q2 - complex<double>(m_2() * m_2(), -m_2() * Gamma_2()));
                }

                // shares the 1^- pole of a_V, as the tensor form factor does
                complex<double>
                a_T(const complex<double> & m) const
                {
                    return std::polar(abs_T_0_over_S_0() * S_0(), arg_T_0()) / (1.0 - m * m / (m_Bsst() * m_Bsst()));
                }

                complex<double>
                a_S(const complex<double> & m) const
                {
                    return S_0() / (1.0 - m * m / (m_Bs0() * m_Bs0()));
                }

                // the coefficient c_j of the singular term of Phi |a_V|^2 at mu_j
                complex<double>
                coefficient_V(const std::size_t & j) const
                {
                    const complex<double> mu    = poles()[j];
                    const complex<double> rho   = residues()[j] / (2.0 * mu);
                    const complex<double> a_reg = a_V(std::conj(mu)) - rho / (std::conj(mu) - mu);

                    return phase_space(mu) * (2.0 * std::conj(a_reg) * rho + 1.0i * std::norm(rho) / (-mu.imag()));
                }

                double
                pdf(const double & m, const double & c) const
                {
                    if ((m < m_min()) || (m > m_max()) || (std::abs(c) > 1.0))
                    {
                        return 0.0;
                    }

                    const complex<double> a_v = a_V(m), a_t = a_T(m), a_s = a_S(m);

                    return phase_space(m).real() * (0.75 * std::norm(a_v) * (1.0 - c * c) + 0.5 * std::norm(a_s + a_t * c));
                }

                double
                normalization(const double & m_lo, const double & m_hi, const double & c_lo, const double & c_hi) const
                {
                    const double lo = std::max(m_lo, m_min()), hi = std::min(m_hi, m_max());
                    const double cl = std::max(c_lo, -1.0), ch = std::min(c_hi, +1.0);
                    if ((lo >= hi) || (cl >= ch))
                    {
                        return 0.0;
                    }

                    // the integrals of the angular terms over [cl, ch]
                    const double I_V = 0.75 * ((ch - cl) - (power_of<3>(ch) - power_of<3>(cl)) / 3.0);
                    const double I_S = 0.5 * (ch - cl);
                    const double I_X = 0.5 * (ch * ch - cl * cl);
                    const double I_T = 0.5 * (power_of<3>(ch) - power_of<3>(cl)) / 3.0;

                    const auto                           mu = poles();
                    const std::array<complex<double>, 2> c{ coefficient_V(0), coefficient_V(1) };

                    const auto remainder = [&](const double & m) -> double
                    {
                        const complex<double> a_v = a_V(m), a_t = a_T(m), a_s = a_S(m);
                        double                result = I_V * (phase_space(m).real() * std::norm(a_v) - (c[0] / (m - mu[0])).real() - (c[1] / (m - mu[1])).real());

                        return result + phase_space(m).real() * (I_S * std::norm(a_s) + I_X * (a_s * std::conj(a_t)).real() + I_T * std::norm(a_t));
                    };

                    double singular = 0.0;
                    for (std::size_t j = 0; j < 2; ++j)
                    {
                        singular += I_V * (c[j] * (std::log(hi - mu[j]) - std::log(lo - mu[j]))).real();
                    }

                    // The remainder is smooth, but near the poles its cancellation leaves noise of about 1e-9 relative,
                    // which a relative tolerance cannot overcome on a short range; the result needs only a tolerance
                    // relative to the singular part, which dominates there.
                    const auto config = GSL::QAGS::Config().epsrel(1.0e-8).epsabs(1.0e-10 * std::abs(singular));

                    return integrate<GSL::QAGS>(remainder, lo, hi, config) + singular;
                }

                double
                re_mu_1() const
                {
                    return poles()[0].real();
                }

                double
                im_mu_1() const
                {
                    return poles()[0].imag();
                }

                double
                re_mu_2() const
                {
                    return poles()[1].real();
                }

                double
                im_mu_2() const
                {
                    return poles()[1].imag();
                }

                // the coefficients of the singular terms of P at mu_j
                double
                re_c_1(const double & c) const
                {
                    return (std::abs(c) > 1.0) ? 0.0 : 0.75 * (1.0 - c * c) * coefficient_V(0).real();
                }

                double
                im_c_1(const double & c) const
                {
                    return (std::abs(c) > 1.0) ? 0.0 : 0.75 * (1.0 - c * c) * coefficient_V(0).imag();
                }

                double
                re_c_2(const double & c) const
                {
                    return (std::abs(c) > 1.0) ? 0.0 : 0.75 * (1.0 - c * c) * coefficient_V(1).real();
                }

                double
                im_c_2(const double & c) const
                {
                    return (std::abs(c) > 1.0) ? 0.0 : 0.75 * (1.0 - c * c) * coefficient_V(1).imag();
                }

                static ObservableEntry::OptionIterator
                begin_options()
                {
                    return options.begin();
                }

                static ObservableEntry::OptionIterator
                end_options()
                {
                    return options.end();
                }
        };

        const std::vector<OptionSpecification> BenchmarkB2::options{};

        const std::set<ReferenceName> BenchmarkB2::references{};

        const std::string BenchmarkB2::description = "Benchmark B2: B -> K mu mu-like distribution in m_mumu and cos(theta_l) with two narrow resonances.";
    } // namespace

    void
    register_benchmarks()
    {
        // a function-local static is initialized exactly once, also when called concurrently
        static const bool registered = []() -> bool
        {
            using std::literals::string_literals::operator""s;

            // benchmark B1
            {
                const std::vector<std::tuple<const char *, const char *, Unit, double, double, double>> parameters{
                    {      "B1::x_1@RvDV:2026A",      R"(x_1)",  Unit::GeV(),               3.1,   3.0,    3.2 },
                    {  "B1::gamma_1@RvDV:2026A", R"(\gamma_1)",  Unit::GeV(),            5.0e-5,   0.0, 1.0e-3 },
                    {      "B1::x_2@RvDV:2026A",      R"(x_2)",  Unit::GeV(),               3.7,   3.6,    3.8 },
                    {  "B1::gamma_2@RvDV:2026A", R"(\gamma_2)",  Unit::GeV(),            1.5e-4,   0.0, 1.0e-3 },
                    { "B1::abs{R_1}@RvDV:2026A",    R"(|R_1|)",  Unit::GeV(),          0.099047,   0.0,    1.0 },
                    { "B1::arg{R_1}@RvDV:2026A", R"(\arg R_1)", Unit::None(),       M_PI / 10.0, -M_PI,   M_PI },
                    { "B1::abs{R_2}@RvDV:2026A",    R"(|R_2|)",  Unit::GeV(),          0.049153,   0.0,    1.0 },
                    { "B1::arg{R_2}@RvDV:2026A", R"(\arg R_2)", Unit::None(), -4.0 * M_PI / 5.0, -M_PI,   M_PI },
                    {        "B1::W@RvDV:2026A",        R"(W)",  Unit::GeV(),               4.6,   4.0,    5.0 },
                };
                for (const auto & [name, latex, unit, central, min, max] : parameters)
                {
                    Parameters::declare(name, latex, unit, central, min, max);
                }

                make_observable("B1::P(m_mumu)@RvDV:2026A", Unit::None(), &BenchmarkB1::pdf, std::make_tuple("m_mumu"));
                make_observable("B1::N@RvDV:2026A", Unit::GeV(), &BenchmarkB1::normalization, std::make_tuple("m_mumu_min", "m_mumu_max"));
                make_observable("B1::Re{mu_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::re_mu_1);
                make_observable("B1::Im{mu_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::im_mu_1);
                make_observable("B1::Re{mu_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::re_mu_2);
                make_observable("B1::Im{mu_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::im_mu_2);
                make_observable("B1::Re{c_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::re_c_1);
                make_observable("B1::Im{c_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::im_c_1);
                make_observable("B1::Re{c_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::re_c_2);
                make_observable("B1::Im{c_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB1::im_c_2);

                auto signal_pdf = make_signal_pdf("benchmark::B1@RvDV:2026A",
                                                  "Benchmark B1: a constant amplitude interfering with two narrow poles within a box in m_mumu.",
                                                  Options{},
                                                  "B1::P(m_mumu)@RvDV:2026A",
                                                  std::make_tuple("m_mumu"s),
                                                  "B1::N@RvDV:2026A",
                                                  std::make_tuple("m_mumu_min"s, "m_mumu_max"s));
                SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
            }

            // benchmark B2
            {
                const std::vector<std::tuple<const char *, const char *, Unit, double, double, double>> parameters{
                    {          "B2::m_B@RvDV:2026A",        R"(m_B)",  Unit::GeV(),             5.279,    5.0,    5.5 },
                    {          "B2::m_K@RvDV:2026A",        R"(m_K)",  Unit::GeV(),             0.494,    0.4,    0.6 },
                    {         "B2::m_mu@RvDV:2026A",      R"(m_\mu)",  Unit::GeV(),             0.106,    0.1,   0.11 },
                    {      "B2::m_B_s^*@RvDV:2026A",  R"(m_{B_s^*})",  Unit::GeV(),             5.415,    5.0,    6.0 },
                    {       "B2::m_B_s0@RvDV:2026A", R"(m_{B_{s0}})",  Unit::GeV(),             5.699,    5.0,    6.5 },
                    {          "B2::m_1@RvDV:2026A",        R"(m_1)",  Unit::GeV(),               3.1,    3.0,    3.2 },
                    {      "B2::Gamma_1@RvDV:2026A",   R"(\Gamma_1)",  Unit::GeV(),            1.0e-4, 1.0e-6, 1.0e-2 },
                    {          "B2::m_2@RvDV:2026A",        R"(m_2)",  Unit::GeV(),               3.7,    3.6,    3.8 },
                    {      "B2::Gamma_2@RvDV:2026A",   R"(\Gamma_2)",  Unit::GeV(),            3.0e-4, 1.0e-6, 1.0e-2 },
                    {          "B2::V_0@RvDV:2026A",        R"(V_0)", Unit::None(),               1.0,    0.0,    2.0 },
                    {          "B2::S_0@RvDV:2026A",        R"(S_0)", Unit::None(),               0.3,    0.0,    1.0 },
                    { "B2::abs{T_0}/S_0@RvDV:2026A",  R"(|T_0|/S_0)", Unit::None(),         1.0 / 3.0,    0.0,    1.0 },
                    {     "B2::arg{T_0}@RvDV:2026A",   R"(\arg T_0)", Unit::None(),        M_PI / 6.0,  -M_PI,   M_PI },
                    {     "B2::abs{R_1}@RvDV:2026A",      R"(|R_1|)", Unit::None(),             0.784,    0.0,    2.0 },
                    {     "B2::arg{R_1}@RvDV:2026A",   R"(\arg R_1)", Unit::None(),       M_PI / 10.0,  -M_PI,   M_PI },
                    {     "B2::abs{R_2}@RvDV:2026A",      R"(|R_2|)", Unit::None(),             0.645,    0.0,    2.0 },
                    {     "B2::arg{R_2}@RvDV:2026A",   R"(\arg R_2)", Unit::None(), -4.0 * M_PI / 5.0,  -M_PI,   M_PI },
                };
                for (const auto & [name, latex, unit, central, min, max] : parameters)
                {
                    Parameters::declare(name, latex, unit, central, min, max);
                }

                make_observable("B2::P(m_mumu,cos(theta_l))@RvDV:2026A", Unit::None(), &BenchmarkB2::pdf, std::make_tuple("m_mumu", "cos(theta_l)"));
                make_observable("B2::N@RvDV:2026A", Unit::None(), &BenchmarkB2::normalization, std::make_tuple("m_mumu_min", "m_mumu_max", "cos(theta_l)_min", "cos(theta_l)_max"));
                make_observable("B2::Re{mu_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB2::re_mu_1);
                make_observable("B2::Im{mu_1}@RvDV:2026A", Unit::GeV(), &BenchmarkB2::im_mu_1);
                make_observable("B2::Re{mu_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB2::re_mu_2);
                make_observable("B2::Im{mu_2}@RvDV:2026A", Unit::GeV(), &BenchmarkB2::im_mu_2);
                make_observable("B2::Re{c_1}(cos(theta_l))@RvDV:2026A", Unit::None(), &BenchmarkB2::re_c_1, std::make_tuple("cos(theta_l)"));
                make_observable("B2::Im{c_1}(cos(theta_l))@RvDV:2026A", Unit::None(), &BenchmarkB2::im_c_1, std::make_tuple("cos(theta_l)"));
                make_observable("B2::Re{c_2}(cos(theta_l))@RvDV:2026A", Unit::None(), &BenchmarkB2::re_c_2, std::make_tuple("cos(theta_l)"));
                make_observable("B2::Im{c_2}(cos(theta_l))@RvDV:2026A", Unit::None(), &BenchmarkB2::im_c_2, std::make_tuple("cos(theta_l)"));

                auto signal_pdf = make_signal_pdf("benchmark::B2@RvDV:2026A",
                                                  "Benchmark B2: B -> K mu mu-like distribution in m_mumu and cos(theta_l) with two narrow resonances.",
                                                  Options{},
                                                  "B2::P(m_mumu,cos(theta_l))@RvDV:2026A",
                                                  std::make_tuple("m_mumu"s, "cos(theta_l)"s),
                                                  "B2::N@RvDV:2026A",
                                                  std::make_tuple("m_mumu_min"s, "m_mumu_max"s, "cos(theta_l)_min"s, "cos(theta_l)_max"s));
                SignalPDFEntries::instance()->insert_or_assign(signal_pdf.first, signal_pdf.second);
            }

            return true;
        }();
        (void) registered;
    }
} // namespace eos
