/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
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
#include <eos/nonleptonic-amplitudes/su3f-amplitudes.hh>

using namespace test;
using namespace eos;

class SU3AmplitudesTest :
    public TestCase
{
    public:
        SU3AmplitudesTest() :
            TestCase("su3_amplitudes_test")
        {
        }

        virtual void run() const
        {

            /* Test with real amplitudes */
            {
                Parameters p = Parameters::Defaults();
                p["nonleptonic::Re{AT3}@SU3F"]  =  0.1;
                p["nonleptonic::Re{CT3}@SU3F"]  = -0.2;
                p["nonleptonic::Re{AT6}@SU3F"]  =  0.3;
                p["nonleptonic::Re{CT6}@SU3F"]  = -0.4;
                p["nonleptonic::Re{AT15}@SU3F"] =  0.5;
                p["nonleptonic::Re{CT15}@SU3F"] = -0.6;
                p["nonleptonic::Re{BT3}@SU3F"]  =  0.7;
                p["nonleptonic::Re{BT6}@SU3F"]  = -0.8;
                p["nonleptonic::Re{BT15}@SU3F"] =  0.9;
                p["nonleptonic::Re{DT3}@SU3F"]  = -1.0;
                p["nonleptonic::Re{AP3}@SU3F"]  =  1.1;
                p["nonleptonic::Re{CP3}@SU3F"]  = -1.2;
                p["nonleptonic::Re{AP6}@SU3F"]  =  1.3;
                p["nonleptonic::Re{CP6}@SU3F"]  = -1.4;
                p["nonleptonic::Re{AP15}@SU3F"] =  1.5;
                p["nonleptonic::Re{CP15}@SU3F"] = -1.6;
                p["nonleptonic::Re{BP3}@SU3F"]  =  1.7;
                p["nonleptonic::Re{BP6}@SU3F"]  = -1.8;
                p["nonleptonic::Re{BP15}@SU3F"] =  1.9;
                p["nonleptonic::Re{DP3}@SU3F"]  = -2.0;
                p["eta::theta_18"]              =  0.0;

                static const double eps = 1.0e-6;

                Options o
                {
                    { "q", "u" },
                    { "P1", "pi^0" },
                    { "P2", "pi^+" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.tree_amplitude(),    complex<double>(-0.003701600698,  0.009833235011), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.penguin_amplitude(), complex<double>(-0.076644441470, -0.030849505280), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>( 1.733325903e-7, -6.626574282e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(-1.929579753e-8,  1.479410207e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(),         complex<double>( 1.540367927e-7, -6.478633261e-7), eps);

                Options oo
                {
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> dd(p, oo);

                TEST_CHECK_RELATIVE_ERROR_C(dd.tree_amplitude(),    complex<double>(-0.001121754409,  0.002979920210), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.penguin_amplitude(), complex<double>(-0.015601822450, -0.006279757473), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.ordered_amplitude(), complex<double>( 2.721554933e-8, -1.379284173e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.inverse_amplitude(), complex<double>( 5.446508563e-8, -2.212538491e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.amplitude(),         complex<double>( 8.168063496e-8, -3.591822664e-7), eps);

                Options ooo
                {
                    { "q", "d" },
                    { "P1", "pi^0" },
                    { "P2", "pi^0" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> ddd(p, ooo);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.tree_amplitude(),    complex<double>( 0.001744951303, -0.004635431438), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.penguin_amplitude(), complex<double>( 0.032024793450,  0.012890028500), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.ordered_amplitude(), complex<double>(-6.808014322e-8,  2.785174180e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.inverse_amplitude(), complex<double>(-6.808014322e-8,  2.785174180e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.amplitude(),         complex<double>(-1.361602864e-7,  5.570348360e-7), eps);

                Options oooo
                {
                    { "q", "s" },
                    { "P1", "eta" },
                    { "P2", "Kbar_d" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> dddd(p, oooo);

                TEST_CHECK_RELATIVE_ERROR_C(dddd.tree_amplitude(),    complex<double>(-0.001221211520,  0.003244126218), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.penguin_amplitude(), complex<double>(-0.021454879470, -0.008635621905), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.ordered_amplitude(), complex<double>( 4.446659188e-8, -1.870220329e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.inverse_amplitude(), complex<double>(-3.336976551e-8,  1.084214168e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.amplitude(),         complex<double>( 1.109682636e-8, -7.860061606e-8), eps);
            }

            /* Test with complex amplitudes */
            {
                Parameters p = Parameters::Defaults();
                p["nonleptonic::Re{AT3}@SU3F"]  =  0.1 * std::cos( 0.1);
                p["nonleptonic::Im{AT3}@SU3F"]  =  0.1 * std::sin( 0.1);
                p["nonleptonic::Re{CT3}@SU3F"]  = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{CT3}@SU3F"]  = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{AT6}@SU3F"]  =  0.3 * std::cos( 0.3);
                p["nonleptonic::Im{AT6}@SU3F"]  =  0.3 * std::sin( 0.3);
                p["nonleptonic::Re{CT6}@SU3F"]  = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{CT6}@SU3F"]  = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{AT15}@SU3F"] =  0.5 * std::cos( 0.5);
                p["nonleptonic::Im{AT15}@SU3F"] =  0.5 * std::sin( 0.5);
                p["nonleptonic::Re{CT15}@SU3F"] = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{CT15}@SU3F"] = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{BT3}@SU3F"]  =  0.7 * std::cos( 0.7);
                p["nonleptonic::Im{BT3}@SU3F"]  =  0.7 * std::sin( 0.7);
                p["nonleptonic::Re{BT6}@SU3F"]  = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{BT6}@SU3F"]  = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{BT15}@SU3F"] =  0.9 * std::cos( 0.9);
                p["nonleptonic::Im{BT15}@SU3F"] =  0.9 * std::sin( 0.9);
                p["nonleptonic::Re{DT3}@SU3F"]  = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{DT3}@SU3F"]  = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{AP3}@SU3F"]  =  1.1 * std::cos( 1.1);
                p["nonleptonic::Im{AP3}@SU3F"]  =  1.1 * std::sin( 1.1);
                p["nonleptonic::Re{CP3}@SU3F"]  = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{CP3}@SU3F"]  = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{AP6}@SU3F"]  =  1.3 * std::cos( 1.3);
                p["nonleptonic::Im{AP6}@SU3F"]  =  1.3 * std::sin( 1.3);
                p["nonleptonic::Re{CP6}@SU3F"]  = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{CP6}@SU3F"]  = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{AP15}@SU3F"] =  1.5 * std::cos( 1.5);
                p["nonleptonic::Im{AP15}@SU3F"] =  1.5 * std::sin( 1.5);
                p["nonleptonic::Re{CP15}@SU3F"] = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{CP15}@SU3F"] = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{BP3}@SU3F"]  =  1.7 * std::cos( 1.7);
                p["nonleptonic::Im{BP3}@SU3F"]  =  1.7 * std::sin( 1.7);
                p["nonleptonic::Re{BP6}@SU3F"]  = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{BP6}@SU3F"]  = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{BP15}@SU3F"] =  1.9 * std::cos( 1.9);
                p["nonleptonic::Im{BP15}@SU3F"] =  1.9 * std::sin( 1.9);
                p["nonleptonic::Re{DP3}@SU3F"]  = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{DP3}@SU3F"]  = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]              =  0.0;

                static const double eps = 1.0e-6;

                Options o
                {
                    { "q", "u" },
                    { "P1", "pi^0" },
                    { "P2", "pi^+" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.tree_amplitude(),    complex<double>(-2.273514188e-3,  8.908726610e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.penguin_amplitude(), complex<double>(-1.059348845e-2,  6.225628113e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(-1.248212396e-7, -1.061211560e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(-5.912919404e-7, -9.905965534e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(),         complex<double>(-7.161131800e-7, -2.051808113e-7), eps);

                Options oo
                {
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> dd(p, oo);

                TEST_CHECK_RELATIVE_ERROR_C(dd.tree_amplitude(),    complex<double>(-2.524130408e-3,  1.991137618e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.penguin_amplitude(), complex<double>( 9.005046361e-3, -1.557506871e-2), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.ordered_amplitude(), complex<double>( 1.120340542e-7,  5.345163223e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.inverse_amplitude(), complex<double>(-8.221418647e-7, -2.769768049e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dd.amplitude(),         complex<double>(-7.101078106e-7, -2.235251727e-7), eps);

                Options ooo
                {
                    { "q", "d" },
                    { "P1", "pi^0" },
                    { "P2", "pi^0" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> ddd(p, ooo);

                TEST_CHECK_RELATIVE_ERROR_C(ddd.tree_amplitude(),    complex<double>(-7.873850931e-4, -4.573258399e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.penguin_amplitude(), complex<double>( 4.827637712e-3, -1.377336708e-2), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.ordered_amplitude(), complex<double>( 1.513145804e-7,  3.332215673e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.inverse_amplitude(), complex<double>( 1.513145804e-7,  3.332215673e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(ddd.amplitude(),         complex<double>( 3.026291608e-7,  6.664431346e-8), eps);

                Options oooo
                {
                    { "q", "s" },
                    { "P1", "eta" },
                    { "P2", "Kbar_d" },
                    { "model", "CKM" }
                };

                SU3FRepresentation<PToPP> dddd(p, oooo);

                TEST_CHECK_RELATIVE_ERROR_C(dddd.tree_amplitude(),    complex<double>(-2.628079667e-4,  3.081147729e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.penguin_amplitude(), complex<double>(-1.199965240e-3,  3.303884344e-3), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.ordered_amplitude(), complex<double>(-5.266082582e-8, -1.206428475e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.inverse_amplitude(), complex<double>(-2.113941987e-7, -2.672089228e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(dddd.amplitude(),         complex<double>(-2.640550246e-7, -3.878517703e-8), eps);
            }
        }

} su3_amplitudes_test;
