/*
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/parameters.hh>

#include <test/test.hh>

#include <cmath>
#include <memory>

using namespace test;
using namespace eos;

//////////////////////////////////////////
/// Generic channel with two Pseudoscalars
//////////////////////////////////////////
template <typename T, T... indices>
auto
_parameter_names(const Parameters & p, std::integer_sequence<T, indices...>) -> std::array<eos::Parameter, sizeof...(indices)>
{
    return std::array<eos::Parameter, sizeof...(indices)>{ { p["test::g0_" + std::to_string(indices + 1) + ""]... } };
}

template <unsigned nchannels_, unsigned nresonances_> struct PPchannel : public KMatrix<nchannels_, nresonances_>::Channel
{
        PPchannel(std::string name, Parameter m1, Parameter m2, unsigned l_orbital, const Parameters & p) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, l_orbital, p["test::q0"], _parameter_names(p, std::make_index_sequence<nresonances_>()))
        {
        }

        double mm = this->_m1 - this->_m2;
        double mp = this->_m1 + this->_m2;

        // sqrt of the Källen factor, defined with an absolute value
        double
        sqlk(const complex<double> & s)
        {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }

        const double          pi = M_PI;
        const complex<double> i  = complex<double>(0.0, 1.0);

        complex<double>
        rho(const complex<double> & s)
        {
            if (real(s) < mp * mp)
            {
                return 0.;
            }
            else
            {
                return sqlk(s) / abs(s) / 16.0 / pi;
            }
        }

        complex<double>
        chew_mandelstam(const complex<double> & S)
        {
            double s = real(S);

            complex<double> result = 0.0;
            if (s < mm * mm)
            {
                result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
                result += mp * sqlk(s) * std::log((mm * mm + mp * mp - 2 * s - 2 * sqlk(s)) / (mp * mp - mm * mm));
                result *= -1.0 / (mp * pi * s);

                return result / 16.0 / pi;
            }
            else if (s < mp * mp)
            {
                result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
                result += 2 * mp * sqlk(s) * std::atan(std::sqrt((s - mm * mm) / (mp * mp - s)));
                result *= -1.0 / (mp * pi * s);

                return result / 16.0 / pi;
            }
            else
            {
                result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
                result += mp * sqlk(s) * std::log((2 * s + 2 * sqlk(s) - mm * mm - mp * mp) / (mp * mp - mm * mm));
                result *= -1.0 / (mp * pi * s);
                result += i * sqlk(s) / s;

                return result / 16.0 / pi;
            }
        }
};

template <unsigned nchannels_, unsigned nresonances_> struct resonance : public KMatrix<nchannels_, nresonances_>::Resonance
{
        resonance(std::string name, Parameter m) :
            KMatrix<nchannels_, nresonances_>::Resonance(name, m)
        {
        }
};

complex<double>
BreitWigner(const complex<double> s, const complex<double> M, const complex<double> Ga)
{
    return M * M * Ga * Ga / (power_of<2>(s - M * M) + M * M * Ga * Ga);
}

/////////
/// TESTS
/////////
class KMatrixTest : public TestCase
{
    public:
        KMatrixTest() :
            TestCase("KMatrix tests")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-4;

            Parameters::declare("test::g0_1", "g_0^1", Unit::None(), 1.0);
            Parameters::declare("test::g0_2", "g_0^2", Unit::None(), 2.2);
            Parameters::declare("test::c_1", "c_1", Unit::None(), 0.0);
            Parameters::declare("test::m_1", "m_1", Unit::GeV(), 20.0);
            Parameters::declare("test::m_2", "m_2", Unit::GeV(), 2.0);
            Parameters::declare("test::m_a", "m_a", Unit::GeV(), 0.7);
            Parameters::declare("test::m_b", "m_b", Unit::GeV(), 0.8);
            Parameters::declare("test::q0", "q0", Unit::GeV(), 0.2);

            Parameters p = Parameters::Defaults();

            // One channel, one resonnance:
            // In the limit where the resonance mass is much larger
            // than the channel masses, one recover a simple Breit-Wigner distribution
            {
                auto res  = std::make_shared<resonance<1, 1>>("res", p["test::m_1"]);
                auto chan = std::make_shared<PPchannel<1, 1>>("chan", p["test::m_a"], p["test::m_b"], 1, p);

                KMatrix<1, 1> simplest_kmatrix({ chan }, { res }, { { p["test::c_1"] } }, "simplest_kmatrix");

                TEST_CHECK_EQUAL(res->_m, 20.0);

                auto qres = 0.5
                            * std::sqrt(std::abs(eos::lambda(power_of<2>((double) res->_m), power_of<2>((double) chan->_m1), power_of<2>((double) chan->_m2))
                                                 / power_of<2>((double) res->_m)));
                auto q300 = 0.5 * std::sqrt(std::abs(eos::lambda(300., power_of<2>((double) chan->_m1), power_of<2>((double) chan->_m2)) / 300.));

                const complex<double> BWfactorat300  = kmatrix_utils::blatt_weisskopf_factor(chan->_l_orbital, q300 / chan->_q0);
                const complex<double> BWfactoratmres = kmatrix_utils::blatt_weisskopf_factor(chan->_l_orbital, qres / chan->_q0);

                // Check that the barrier factors are 1 at large q^2 = 300 GeV^2
                TEST_CHECK_NEARLY_EQUAL(std::abs(pow(q300 / qres, chan->_l_orbital) * BWfactorat300 / BWfactoratmres), 1., eps);

                // Check |T-Matrix|^2 against Breit-Wigner
                // The mass of the BW gets correction from the channels loop
                const double          m = res->_m;
                const complex<double> M = m * (1.0 - 0.5 * chan->chew_mandelstam(m * m) / m / m);

                // s = 100. (everybody is ~zero there)
                double Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 100.)[0]);
                TEST_CHECK_NEARLY_EQUAL(std::abs(BreitWigner(100., M, 1. / M)), Trowsq, eps);

                // s = 380. (resonance region)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 380.)[0]);
                TEST_CHECK_NEARLY_EQUAL(std::abs(BreitWigner(380., M, 1. / M)), Trowsq, eps);

                // s = 420. (resonance region)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 420.)[0]);
                TEST_CHECK_NEARLY_EQUAL(std::abs(BreitWigner(420., M, 1. / M)), Trowsq, eps);

                // s = 800. (everybody is ~zero there)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 800.)[0]);
                TEST_CHECK_NEARLY_EQUAL(std::abs(BreitWigner(800., M, 1. / M)), Trowsq, eps);
            }

            p["test::g0_1"] = 1.1;
            p["test::m_1"]  = 1.0;

            // One channel, two resonances
            {
                auto res1 = std::make_shared<resonance<1, 2>>("res1", p["test::m_1"]);
                auto res2 = std::make_shared<resonance<1, 2>>("res2", p["test::m_2"]);

                auto chan_12 = std::make_shared<PPchannel<1, 2>>("chan_12", p["test::m_a"], p["test::m_b"], 1, p);

                // Test the phase space and Chew-Mandelstam functions for a PP channel
                TEST_CHECK_NEARLY_EQUAL(std::abs(chan_12->rho(9.0)), 0.0172195, eps);
                TEST_CHECK_NEARLY_EQUAL(std::abs(chan_12->rho(1.5)), 0., eps);

                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(9.0).real(), -0.0144157, eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(9.0).imag(), 0.0172195, eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(1.5).real(), -0.00854099, eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(1.5).imag(), 0., eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(0.001).real(), -0.0126445, eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->chew_mandelstam(0.001).imag(), 0., eps);

                // Test K matrix inversion into T matrix
                KMatrix<1, 2> kmatrix_12({ chan_12 }, { res1, res2 }, { { p["test::c_1"] } }, "kmatrix_12");

                auto kmatrix_12_ats0 = kmatrix_12.tmatrix_row(0, 9.0);
                auto kmatrix_12_ats1 = kmatrix_12.tmatrix_row(0, 1.5);
                // At the resonance mass
                auto kmatrix_12_ats2 = kmatrix_12.tmatrix_row(0, 1.0);
                auto kmatrix_12_ats3 = kmatrix_12.tmatrix_row(0, 0.001);

                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats0[0].real(), -1.11080, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats0[0].imag(), 0.0217595, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats1[0].real(), -0.618935, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats1[0].imag(), 0, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats2[0].real(), 111.325, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats2[0].imag(), 0, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats3[0].real(), 2.33115, eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats3[0].imag(), 0, eps);

                TEST_CHECK_NEARLY_EQUAL(kmatrix_12.width(1), 0.0318047, eps);
            }
        }
} kmatrix_test;
