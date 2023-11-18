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

#include <test/test.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/parameters.hh>

#include <cmath>
#include <memory>

using namespace test;
using namespace eos;

///////////////////////
/// Simplest K matrix, 1 channel, 1 resonance
///////////////////////
struct PPchan11 :
    public KMatrix<1,1>::Channel
{
    PPchan11(std::string name, double m1, double m2, unsigned N_orbital, const Parameters & p) :
        Channel(name, m1, m2, N_orbital, {{ p["test::g0_1"] }})
    {
    };

    // Useful definitions for beta and rho
    double mm = this->_m1 - this->_m2;
    double mp = this->_m1 + this->_m2;
    // sqrt of the Källen factor, defined with an absolute value
    double sqlk(const double & s)
    {
        return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
    }
    const double pi = M_PI;
    const complex<double> i = complex<double>(0.0, 1.0);

    double beta(const double & s)
    {
        if (s < mp * mp)
            return 0.;
        else
            return sqlk(s) / s;
    }

    complex<double> rho(const double & s)
    {
        complex<double> result = 0.0;
        if (s < mm * mm)
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += mp * sqlk(s) * std::log((mm * mm + mp * mp - 2 * s - 2 * sqlk(s)) / (mp * mp - mm * mm));
            result *= i / (mp * pi * s);

            return result;
        }
        else if (s < mp * mp)
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += 2 * mp * sqlk(s) * std::atan(std::sqrt((s - mm * mm) / (mp * mp - s)));
            result *= i / (mp * pi * s);

            return result;
        }
        else
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += mp * sqlk(s) * std::log((2 * s + 2 * sqlk(s) - mm * mm - mp * mp) / (mp * mp - mm * mm));
            result *= i / (mp * pi * s);
            result += sqlk(s) / s;

            return result;
        }
    }

};

struct res11 :
    public KMatrix<1,1>::Resonance
{
    res11(std::string name, Parameter m) :
        Resonance(name, m)
    {
    };
};

double BreitWigner(double const s, double const M, double const Ga)
{
    return M * M * Ga * Ga / (pow(s - M * M, 2) + M * M * Ga * Ga);
}


///////////////////////
/// 1 channel, 2 resonances
///////////////////////
struct PPchan12 :
    public KMatrix<1,2>::Channel
{

    PPchan12(std::string name, double m1, double m2, unsigned N_orbital, const Parameters & p) :
        Channel(name, m1, m2, N_orbital, {{ p["test::g0_1"], p["test::g0_2"] }})
    {
    };

    // Useful definitions for beta and rho
    double mm = this->_m1 - this->_m2;
    double mp = this->_m1 + this->_m2;
    // sqrt of the Källen factor, defined with an absolute value
    double sqlk(const double & s)
    {
        return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
    }
    const double pi = M_PI;
    const complex<double> i = complex<double>(0.0, 1.0);

    double beta(const double & s)
    {
        if (s < mp * mp)
            return 0.;
        else
            return sqlk(s) / s;
    }

    complex<double> rho(const double & s)
    {
        complex<double> result = 0.0;
        if (s < mm * mm)
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += mp * sqlk(s) * std::log((mm * mm + mp * mp - 2 * s - 2 * sqlk(s)) / (mp * mp - mm * mm));
            result *= i / (mp * pi * s);

            return result;
        }
        else if (s < mp * mp)
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += 2 * mp * sqlk(s) * std::atan(std::sqrt((s - mm * mm) / (mp * mp - s)));
            result *= i / (mp * pi * s);

            return result;
        }
        else
        {
            result += mm * (mp * mp - s) * std::log((mm + mp) / (mp - mm));
            result += mp * sqlk(s) * std::log((2 * s + 2 * sqlk(s) - mm * mm - mp * mp)/(mp * mp - mm * mm));
            result *= i / (mp * pi * s);
            result += sqlk(s) / s;

            return result;
        }
    }
};

struct res12 :
    public KMatrix<1,2>::Resonance
{
    res12(std::string name, Parameter m) :
        Resonance(name, m)
    {
    };
};


///////////////////////
/// TESTS
///////////////////////

class KMatrixTest :
    public TestCase
{
    public:
        KMatrixTest() :
            TestCase("KMatrix tests")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-4;

            Parameters::declare("test::g0_1", "g_0^1", Unit::None(),  1.0);
            Parameters::declare("test::c_1",  "c_1",   Unit::None(),  0.0);
            Parameters::declare("test::m_1",  "m_1",   Unit::GeV(),  15.0);
            Parameters::declare("test::g0_2", "g_0^2", Unit::None(),  2.2);
            Parameters::declare("test::m_2",  "m_2",   Unit::GeV(),   2.0);

            Parameters p = Parameters::Defaults();

            // One channel, one resonnance:
            // In the limit where the resonance mass is much larger
            // than the channel masses, one recover a simple Breit-Wigner distribution
            {
                auto res = std::make_shared<res11>("res", p["test::m_1"]);
                auto chan = std::make_shared<PPchan11>("chan", 0.7, 0.8, 3, p);

                KMatrix<1,1> simplest_kmatrix({chan}, {res}, {{p["test::c_1"]}}, "simplest_kmatrix");

                // Check |T-Matrix|^2 against Breit-Wigner
                // The mass of the BW gets correction from the channels loop

                double m = res->_m;
                double M = std::sqrt(m*m + chan->rho(m*m).imag());

                // s = 100. (everybody is ~zero there)
                double Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 100.)[0]);
                TEST_CHECK_NEARLY_EQUAL(BreitWigner(100., M, 1./M),        Trowsq,        eps);

                // s = 222. (resonance region)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 222.)[0]);
                TEST_CHECK_NEARLY_EQUAL(BreitWigner(222., M, 1./M),        Trowsq,        eps);

                // s = 235. (resonance region)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 235.)[0]);
                TEST_CHECK_NEARLY_EQUAL(BreitWigner(235., M, 1./M),        Trowsq,        eps);

                // s = 300. (everybody is ~zero there)
                Trowsq = std::norm(simplest_kmatrix.tmatrix_row(0, 300.)[0]);
                TEST_CHECK_NEARLY_EQUAL(BreitWigner(300., M, 1./M),        Trowsq,        eps);
            }

            p["test::g0_1"] = 1.1;
            p["test::m_1"] = 1.0;

            // One channel, two resonances
            {
                auto res1 = std::make_shared<res12>("res1", p["test::m_1"]);
                auto res2 = std::make_shared<res12>("res2", p["test::m_2"]);

                auto chan_12 = std::make_shared<PPchan12>("chan_12", 0.7, 0.8, 3, p);

                // Test beta and rho for a PP channel
                TEST_CHECK_NEARLY_EQUAL(chan_12->beta(9.0),          0.865544,  eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->beta(1.5),          0.,        eps);

                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(9.0).real(),    0.865544,  eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(9.0).imag(),    0.724611,  eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(1.5).real(),    0.,        eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(1.5).imag(),    0.429317,  eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(0.001).real(),  0.,        eps);
                TEST_CHECK_NEARLY_EQUAL(chan_12->rho(0.001).imag(),  0.635582,  eps);

                // Test KMatrix inversion into TMatrix
                KMatrix<1,2> kmatrix_12({chan_12}, {res1, res2}, {{p["test::c_1"]}}, "kmatrix_12");

                auto kmatrix_12_ats0 = kmatrix_12.tmatrix_row(0, 9.0);
                auto kmatrix_12_ats1 = kmatrix_12.tmatrix_row(0, 1.5);
                 // At the resonance mass
                auto kmatrix_12_ats2 = kmatrix_12.tmatrix_row(0, 1.0);
                auto kmatrix_12_ats3 = kmatrix_12.tmatrix_row(0, 0.001);

                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats0[0].real(),  -0.217114,   eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats0[0].imag(),   1.11299,    eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats1[0].real(),  -0.610949,   eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats1[0].imag(),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats2[0].real(),   1.928403,   eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats2[0].imag(),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats3[0].real(),   0.953701,   eps);
                TEST_CHECK_NEARLY_EQUAL(kmatrix_12_ats3[0].imag(),   0.,         eps);

                TEST_CHECK_NEARLY_EQUAL(kmatrix_12.width(1),        1.598678,   eps);
            }
    }
} kmatrix_test;
