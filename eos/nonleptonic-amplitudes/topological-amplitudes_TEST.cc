/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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
#include <eos/nonleptonic-amplitudes/topological-amplitudes.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class TopologicalAmplitudesTest : public TestCase
{
    public:
        TopologicalAmplitudesTest() :
            TestCase("topological_amplitudes_test")
        {
        }

        virtual void
        run() const
        {
            /* Test with real amplitudes */
            {
                Parameters p                          = Parameters::Defaults();
                p["nonleptonic::Re{T}@Topological"]   = 0.1;
                p["nonleptonic::Re{C}@Topological"]   = -0.2;
                p["nonleptonic::Re{A}@Topological"]   = 0.3;
                p["nonleptonic::Re{E}@Topological"]   = -0.4;
                p["nonleptonic::Re{TES}@Topological"] = 0.5;
                p["nonleptonic::Re{TAS}@Topological"] = -0.6;
                p["nonleptonic::Re{TS}@Topological"]  = 0.7;
                p["nonleptonic::Re{TPA}@Topological"] = -0.8;
                p["nonleptonic::Re{TP}@Topological"]  = 0.9;
                p["nonleptonic::Re{TSS}@Topological"] = -1.0;
                p["nonleptonic::Re{P}@Topological"]   = 1.1;
                p["nonleptonic::Re{PT}@Topological"]  = -1.2;
                p["nonleptonic::Re{S}@Topological"]   = 1.3;
                p["nonleptonic::Re{PC}@Topological"]  = -1.4;
                p["nonleptonic::Re{PTA}@Topological"] = 1.5;
                p["nonleptonic::Re{PA}@Topological"]  = -1.6;
                p["nonleptonic::Re{PTE}@Topological"] = 1.7;
                p["nonleptonic::Re{PAS}@Topological"] = -1.8;
                p["nonleptonic::Re{PSS}@Topological"] = 1.9;
                p["nonleptonic::Re{PES}@Topological"] = -2.0;
                p["eta::theta_18"]                    = 0.0;

                static const double eps = 1.0e-6;

                Options o{
                    {            "q",    "u" },
                    {           "P1", "pi^0" },
                    {           "P2", "pi^+" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.tree_amplitude(), complex<double>(0.001145733549, -0.003043620361), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.penguin_amplitude(), complex<double>(0.008128955914, 0.003271917227), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(-1.882888193e-9, 7.649339906e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(5.006745077e-8, -2.017304058e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(), complex<double>(4.818456257e-8, -1.252370068e-7), eps);

                Options oo{
                    {            "q",    "d" },
                    {           "P1", "pi^+" },
                    {           "P2", "pi^-" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> dd(p, oo);

                TEST_CHECK_RELATIVE_ERROR_C(dd.tree_amplitude(), complex<double>(-0.001495672546, 0.003973226947), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.penguin_amplitude(), complex<double>(0.000821148550, 0.000330513551), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.ordered_amplitude(), complex<double>(-3.549528431e-8, -5.563165578e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.inverse_amplitude(), complex<double>(5.180227961e-8, -1.130758467e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.amplitude(), complex<double>(1.630699530e-8, -1.186390123e-7), eps);

                Options ooo{
                    {            "q",    "d" },
                    {           "P1", "pi^0" },
                    {           "P2", "pi^0" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> ddd(p, ooo);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.tree_amplitude(), complex<double>(-0.000560877204, 0.001489960105), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.penguin_amplitude(), complex<double>(0.004105742750, 0.001652567756), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.ordered_amplitude(), complex<double>(-2.591813329e-8, 2.923643060e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.inverse_amplitude(), complex<double>(-2.591813329e-8, 2.923643060e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.amplitude(), complex<double>(-5.183626659e-8, 5.847286120e-8), eps);

                Options oooo{
                    {            "q",      "s" },
                    {           "P1",    "eta" },
                    {           "P2", "Kbar_d" },
                    {        "model",    "CKM" },
                    { "cp-conjugate",   "true" }
                };

                TopologicalRepresentation<PToPP> dddd(p, oooo);

                TEST_CHECK_RELATIVE_ERROR_C(dddd.tree_amplitude(), complex<double>(-0.000915908639, 0.002433094663), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.penguin_amplitude(), complex<double>(-0.007375114819, -0.002968495030), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.ordered_amplitude(), complex<double>(4.415737482e-9, -6.838057151e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.inverse_amplitude(), complex<double>(1.114241606e-8, -5.356868028e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.amplitude(), complex<double>(1.555815354e-8, -7.373743954e-8), eps);
            }

            /* Test with complex amplitudes and theta = 0*/
            {
                Parameters p                          = Parameters::Defaults();
                p["nonleptonic::Re{T}@Topological"]   = 0.1 * std::cos(0.1);
                p["nonleptonic::Im{T}@Topological"]   = 0.1 * std::sin(0.1);
                p["nonleptonic::Re{C}@Topological"]   = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{C}@Topological"]   = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{A}@Topological"]   = 0.3 * std::cos(0.3);
                p["nonleptonic::Im{A}@Topological"]   = 0.3 * std::sin(0.3);
                p["nonleptonic::Re{E}@Topological"]   = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{E}@Topological"]   = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{TES}@Topological"] = 0.5 * std::cos(0.5);
                p["nonleptonic::Im{TES}@Topological"] = 0.5 * std::sin(0.5);
                p["nonleptonic::Re{TAS}@Topological"] = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{TAS}@Topological"] = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{TS}@Topological"]  = 0.7 * std::cos(0.7);
                p["nonleptonic::Im{TS}@Topological"]  = 0.7 * std::sin(0.7);
                p["nonleptonic::Re{TPA}@Topological"] = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{TPA}@Topological"] = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{TP}@Topological"]  = 0.9 * std::cos(0.9);
                p["nonleptonic::Im{TP}@Topological"]  = 0.9 * std::sin(0.9);
                p["nonleptonic::Re{TSS}@Topological"] = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{TSS}@Topological"] = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{P}@Topological"]   = 1.1 * std::cos(1.1);
                p["nonleptonic::Im{P}@Topological"]   = 1.1 * std::sin(1.1);
                p["nonleptonic::Re{PT}@Topological"]  = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{PT}@Topological"]  = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{S}@Topological"]   = 1.3 * std::cos(1.3);
                p["nonleptonic::Im{S}@Topological"]   = 1.3 * std::sin(1.3);
                p["nonleptonic::Re{PC}@Topological"]  = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{PC}@Topological"]  = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{PTA}@Topological"] = 1.5 * std::cos(1.5);
                p["nonleptonic::Im{PTA}@Topological"] = 1.5 * std::sin(1.5);
                p["nonleptonic::Re{PA}@Topological"]  = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{PA}@Topological"]  = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{PTE}@Topological"] = 1.7 * std::cos(1.7);
                p["nonleptonic::Im{PTE}@Topological"] = 1.7 * std::sin(1.7);
                p["nonleptonic::Re{PAS}@Topological"] = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{PAS}@Topological"] = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{PSS}@Topological"] = 1.9 * std::cos(1.9);
                p["nonleptonic::Im{PSS}@Topological"] = 1.9 * std::sin(1.9);
                p["nonleptonic::Re{PES}@Topological"] = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{PES}@Topological"] = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                    = 0.0;

                static const double eps = 1.0e-6;

                Options o{
                    {            "q",    "u" },
                    {           "P1", "pi^0" },
                    {           "P2", "pi^+" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.tree_amplitude(), complex<double>(2.714849536e-3, -1.505497114e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.penguin_amplitude(), complex<double>(-7.413420944e-3, 2.127194598e-2), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(-1.630246346e-7, -3.875166918e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(5.413850131e-8, -4.135936105e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(), complex<double>(-1.088861333e-7, -8.011103022e-8), eps);

                Options oo{
                    {            "q",    "d" },
                    {           "P1", "pi^+" },
                    {           "P2", "pi^-" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> dd(p, oo);

                TEST_CHECK_RELATIVE_ERROR_C(dd.tree_amplitude(), complex<double>(1.261996239e-3, 3.974744683e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.penguin_amplitude(), complex<double>(-1.227283279e-2, 2.640641563e-2), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.ordered_amplitude(), complex<double>(-2.505699225e-7, -9.081234659e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.inverse_amplitude(), complex<double>(-2.639484094e-7, -5.705999870e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.amplitude(), complex<double>(-5.145183319e-7, -1.478723453e-7), eps);

                Options ooo{
                    {            "q",    "d" },
                    {           "P1", "pi^0" },
                    {           "P2", "pi^0" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> ddd(p, ooo);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.tree_amplitude(), complex<double>(2.805870338e-3, 2.431652123e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.penguin_amplitude(), complex<double>(-4.902145577e-3, 1.942516576e-2), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.ordered_amplitude(), complex<double>(-1.802650427e-7, -1.728911992e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.inverse_amplitude(), complex<double>(-1.802650427e-7, -1.728911992e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.amplitude(), complex<double>(-3.605300854e-7, -3.457823985e-8), eps);

                Options oooo{
                    {            "q",      "s" },
                    {           "P1",    "eta" },
                    {           "P2", "Kbar_d" },
                    {        "model",    "CKM" },
                    { "cp-conjugate",   "true" }
                };

                TopologicalRepresentation<PToPP> dddd(p, oooo);

                TEST_CHECK_RELATIVE_ERROR_C(dddd.tree_amplitude(), complex<double>(-2.475246464e-3, 7.949800121e-4), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.penguin_amplitude(), complex<double>(-6.997788583e-4, -7.919254440e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.ordered_amplitude(), complex<double>(5.875775885e-8, -2.618615750e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.inverse_amplitude(), complex<double>(-6.722738747e-8, -9.218878246e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.amplitude(), complex<double>(-8.469628623e-9, -3.540503574e-8), eps);
            }

            /* Test with complex amplitudes and theta = 1.0 */
            {
                Parameters p                          = Parameters::Defaults();
                p["nonleptonic::Re{T}@Topological"]   = 0.1 * std::cos(0.1);
                p["nonleptonic::Im{T}@Topological"]   = 0.1 * std::sin(0.1);
                p["nonleptonic::Re{C}@Topological"]   = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{C}@Topological"]   = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{A}@Topological"]   = 0.3 * std::cos(0.3);
                p["nonleptonic::Im{A}@Topological"]   = 0.3 * std::sin(0.3);
                p["nonleptonic::Re{E}@Topological"]   = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{E}@Topological"]   = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{TES}@Topological"] = 0.5 * std::cos(0.5);
                p["nonleptonic::Im{TES}@Topological"] = 0.5 * std::sin(0.5);
                p["nonleptonic::Re{TAS}@Topological"] = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{TAS}@Topological"] = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{TS}@Topological"]  = 0.7 * std::cos(0.7);
                p["nonleptonic::Im{TS}@Topological"]  = 0.7 * std::sin(0.7);
                p["nonleptonic::Re{TPA}@Topological"] = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{TPA}@Topological"] = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{TP}@Topological"]  = 0.9 * std::cos(0.9);
                p["nonleptonic::Im{TP}@Topological"]  = 0.9 * std::sin(0.9);
                p["nonleptonic::Re{TSS}@Topological"] = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{TSS}@Topological"] = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{P}@Topological"]   = 1.1 * std::cos(1.1);
                p["nonleptonic::Im{P}@Topological"]   = 1.1 * std::sin(1.1);
                p["nonleptonic::Re{PT}@Topological"]  = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{PT}@Topological"]  = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{S}@Topological"]   = 1.3 * std::cos(1.3);
                p["nonleptonic::Im{S}@Topological"]   = 1.3 * std::sin(1.3);
                p["nonleptonic::Re{PC}@Topological"]  = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{PC}@Topological"]  = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{PTA}@Topological"] = 1.5 * std::cos(1.5);
                p["nonleptonic::Im{PTA}@Topological"] = 1.5 * std::sin(1.5);
                p["nonleptonic::Re{PA}@Topological"]  = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{PA}@Topological"]  = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{PTE}@Topological"] = 1.7 * std::cos(1.7);
                p["nonleptonic::Im{PTE}@Topological"] = 1.7 * std::sin(1.7);
                p["nonleptonic::Re{PAS}@Topological"] = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{PAS}@Topological"] = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{PSS}@Topological"] = 1.9 * std::cos(1.9);
                p["nonleptonic::Im{PSS}@Topological"] = 1.9 * std::sin(1.9);
                p["nonleptonic::Re{PES}@Topological"] = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{PES}@Topological"] = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                    = 1.0;

                static const double eps = 1.0e-6;

                Options o{
                    {            "q",    "u" },
                    {           "P1",  "eta" },
                    {           "P2", "pi^+" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.tree_amplitude(), complex<double>(-0.001018379427023484, 0.0005647337976274578), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.penguin_amplitude(), complex<double>(0.0027808816931099217, -0.007979415385035703), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(6.115290434864023e-8, 1.4536313019197382e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(3.9109708177886845e-7, -2.3065477505584753e-8), eps);

                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(), complex<double>(4.5224998612750857e-7, -8.529164486387376e-9), eps);

                Options oo{
                    {            "q",         "d" },
                    {           "P1", "eta_prime" },
                    {           "P2",       "K_d" },
                    {        "model",       "CKM" },
                    { "cp-conjugate",      "true" }
                };

                TopologicalRepresentation<PToPP> dd(p, oo);

                TEST_CHECK_RELATIVE_ERROR_C(dd.tree_amplitude(), complex<double>(0.00045961421670513524, -0.00014761524595177606), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.penguin_amplitude(), complex<double>(-0.013955444602458823, -0.026203892903736634), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.ordered_amplitude(), complex<double>(2.1733519347033603e-7, -1.1130743983491132e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.inverse_amplitude(), complex<double>(5.80912303675578e-7, -3.216600683041854e-9), eps);

                TEST_CHECK_RELATIVE_ERROR_C(dd.amplitude(), complex<double>(7.98247497145914e-7, -1.1452404051795318e-7), eps);


                Options ooo{
                    {            "q",    "d" },
                    {           "P1",  "eta" },
                    {           "P2",  "eta" },
                    {        "model",  "CKM" },
                    { "cp-conjugate", "true" }
                };

                TopologicalRepresentation<PToPP> ddd(p, ooo);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.tree_amplitude(), complex<double>(0.006932851020723416, 0.008997971640539072), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.penguin_amplitude(), complex<double>(-0.03059460579418031, 0.053191625860970544), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.ordered_amplitude(), complex<double>(-5.129113722089845e-7, -1.9515133715782058e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.inverse_amplitude(), complex<double>(-5.129113722089845e-7, -1.9515133715782058e-7), eps);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.amplitude(), complex<double>(-1.0258227444179687e-6, -3.90302674315641e-7), eps);
            }
        }
} topological_amplitudes_test;
