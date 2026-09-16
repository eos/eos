/*
 * Copyright (c) 2026    Carolina Bolognani
 * Copyright (c) 2026    Dominik Suelmann
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

#include <eos/rare-c-decays/d-to-psd-l-l.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class DToPseudoscalarLeptonLeptonBGHT2019Test : public TestCase
{
    public:
        DToPseudoscalarLeptonLeptonBGHT2019Test() :
            TestCase("d_to_psdoscalar_lepton_lepton_BGHT2019_test")
        {
        }

        virtual void
        run() const
        {
            Parameters p                   = Parameters::Defaults();
            p["D->pi::alpha^f+_0@BSZ2015"] = 0.612;
            p["D->pi::alpha^f+_1@BSZ2015"] = -0.9331;
            p["D->pi::alpha^f+_2@BSZ2015"] = -3.3364;

            p["D->pi::alpha^f0_1@BSZ2015"] = 0.2498;
            p["D->pi::alpha^f0_2@BSZ2015"] = -0.2269;

            p["D->pi::alpha^fT_0@BSZ2015"] = 0.5315;
            p["D->pi::alpha^fT_1@BSZ2015"] = 0.0743;
            p["D->pi::alpha^fT_2@BSZ2015"] = -0.9102;

            p["D_d->pi::res_a_rho@BGHT2019"]             = +0.18;
            p["D_d->pi::res_a_omega@BGHT2019"]           = +0.06;
            p["D_d->pi::res_a_phi@BGHT2019"]             = +0.23;
            p["D_d->pi::res_delta_rho@BGHT2019"]         = +0.00;
            p["D_d->pi::res_delta_omega_m_rho@BGHT2019"] = +M_PI;
            p["D_d->pi::res_delta_phi_m_rho@BGHT2019"]   = +0.00;

            p["D_d->pi::res_a_eta@BGHT2019"]                 = +0.00057;
            p["D_d->pi::res_a_eta_prime@BGHT2019"]           = +0.0008;
            p["D_d->pi::res_delta_eta@BGHT2019"]             = +M_PI;
            p["D_d->pi::res_delta_eta_prime_m_eta@BGHT2019"] = +M_PI;

            p["QED::alpha_e(m_c)"]    = 7.606410e-03;
            p["WET::G_Fermi"]         = 1.16637880000000e-05;
            p["mass::rho^0"]          = 0.77526;
            p["mass::omega"]          = 0.78266;
            p["mass::phi"]            = 1.019461;
            p["mass::eta"]            = 0.547862;
            p["mass::eta_prime"]      = 0.95778;
            p["mass::mu"]             = 0.1056583755;
            p["life_time::rho^0"]     = 4.465481194029851e-24;
            p["life_time::omega"]     = 7.583086728110599e-23;
            p["life_time::phi"]       = 1.5490984419863498e-22;
            p["life_time::eta"]       = 5.02451853435115e-19;
            p["life_time::eta_prime"] = 3.50112727659575e-21;

            p["mass::D_u"]      = 1.86483;
            p["mass::D_d"]      = 1.86966;
            p["mass::D_s"]      = 1.86483;
            p["mass::K_u"]      = 0.493677;
            p["mass::pi^0"]     = 0.1349768;
            p["mass::pi^+"]     = 0.13957039;
            p["life_time::D_u"] = 4.1035656359102247e-13;
            p["life_time::D_d"] = 1.0332997299843015e-12;
            p["life_time::D_s"] = 5.043769563218392e-13;

            p["uc::Re{c7}"]     = 0.0;
            p["uc::Im{c7}"]     = 0.0;
            p["ucmumu::Re{c9}"] = 0.0;
            p["ucmumu::Im{c9}"] = 0.0;

            Options oo{
                {                "model"_ok,      "WET"_ov },
                {                  "tag"_ok, "BGHT2019"_ov },
                { "nonlocal-formfactors"_ok, "BGHT2019"_ov },
                {         "form-factors"_ok,  "BSZ2015"_ov },
                {                    "l"_ok,       "mu"_ov },
                {                    "q"_ok,        "d"_ov }
            };

            static const double eps = 1e-5;
            static const double q2  = 0.9;

            DToPseudoscalarLeptonLepton c(p, oo);

            auto x = c.differential_branching_ratio(q2);

            TEST_CHECK_RELATIVE_ERROR(c.differential_branching_ratio(q2), 8.894089645900448e-08, eps);
            TEST_CHECK_RELATIVE_ERROR(c.double_differential_decay_width(q2, 0), 4.144414237386601e-20, eps);
            TEST_CHECK_NEARLY_EQUAL(c.double_differential_decay_width(q2, 1), 0., 1e-20);

            const complex<double> ckm_factor_ds(-5.383319178425827e-05, 1.427709577227314e-04);

            p["uc::Re{c7}"]  = real(complex<double>(2.2, -1.1) / ckm_factor_ds);
            p["uc::Im{c7}"]  = imag(complex<double>(2.2, -1.1) / ckm_factor_ds);
            p["uc::Re{c7'}"] = real(complex<double>(1.2, 0.5) / ckm_factor_ds);
            p["uc::Im{c7'}"] = imag(complex<double>(1.2, 0.5) / ckm_factor_ds);

            p["ucmumu::Re{c9}"]  = real(complex<double>(0.3, -1.0) / ckm_factor_ds);
            p["ucmumu::Im{c9}"]  = imag(complex<double>(0.3, -1.0) / ckm_factor_ds);
            p["ucmumu::Re{c9'}"] = real(complex<double>(-0.6, 0.2) / ckm_factor_ds);
            p["ucmumu::Im{c9'}"] = imag(complex<double>(-0.6, 0.2) / ckm_factor_ds);

            p["ucmumu::Re{c10}"]  = real(complex<double>(2.0, 1.4) / ckm_factor_ds);
            p["ucmumu::Im{c10}"]  = imag(complex<double>(2.0, 1.4) / ckm_factor_ds);
            p["ucmumu::Re{c10'}"] = real(complex<double>(-1.2, 0.4) / ckm_factor_ds);
            p["ucmumu::Im{c10'}"] = imag(complex<double>(-1.2, 0.4) / ckm_factor_ds);
        }
} d_to_psdoscalar_lepton_lepton_BGHT2019_test;
