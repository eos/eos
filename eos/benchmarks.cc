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

            return true;
        }();
        (void) registered;
    }
} // namespace eos
