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

#include <eos/maths/complex.hh>
#include <eos/nonleptonic-amplitudes/qcdf-amplitudes.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class QCDFAmplitudesTest : public TestCase
{
    public:
        QCDFAmplitudesTest() :
            TestCase("qcdf_amplitudes_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p = Parameters::Defaults();

                p["nonleptonic::Re{alpha1}@QCDF"] = 0.1 * std::cos(0.1);
                p["nonleptonic::Im{alpha1}@QCDF"] = 0.1 * std::sin(0.1);
                p["nonleptonic::Re{alpha2}@QCDF"] = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{alpha2}@QCDF"] = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{b2}@QCDF"]     = 0.3 * std::cos(0.3);
                p["nonleptonic::Im{b2}@QCDF"]     = 0.3 * std::sin(0.3);
                p["nonleptonic::Re{b1}@QCDF"]     = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{b1}@QCDF"]     = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{bS2}@QCDF"]    = 0.5 * std::cos(0.5);
                p["nonleptonic::Im{bS2}@QCDF"]    = 0.5 * std::sin(0.5);
                p["nonleptonic::Re{bS1}@QCDF"]    = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{bS1}@QCDF"]    = -0.6 * std::sin(-0.6);

                p["nonleptonic::Re{alpha4_u}@QCDF"] = 0.7 * std::cos(0.7);
                p["nonleptonic::Im{alpha4_u}@QCDF"] = 0.7 * std::sin(0.7);
                p["nonleptonic::Re{alpha3_u}@QCDF"] = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{alpha3_u}@QCDF"] = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{b4_u}@QCDF"]     = 0.9 * std::cos(0.9);
                p["nonleptonic::Im{b4_u}@QCDF"]     = 0.9 * std::sin(0.9);
                p["nonleptonic::Re{bS4_u}@QCDF"]    = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{bS4_u}@QCDF"]    = -1.0 * std::sin(-1.0);

                p["nonleptonic::Re{alpha4EW_c}@QCDF"] = 1.1 * std::cos(1.1);
                p["nonleptonic::Im{alpha4EW_c}@QCDF"] = 1.1 * std::sin(1.1);
                p["nonleptonic::Re{alpha3EW_c}@QCDF"] = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{alpha3EW_c}@QCDF"] = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{b3EW_c}@QCDF"]     = 1.3 * std::cos(1.3);
                p["nonleptonic::Im{b3EW_c}@QCDF"]     = 1.3 * std::sin(1.3);
                p["nonleptonic::Re{b4EW_c}@QCDF"]     = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{b4EW_c}@QCDF"]     = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{bS3EW_c}@QCDF"]    = 1.5 * std::cos(1.5);
                p["nonleptonic::Im{bS3EW_c}@QCDF"]    = 1.5 * std::sin(1.5);
                p["nonleptonic::Re{bS4EW_c}@QCDF"]    = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{bS4EW_c}@QCDF"]    = -1.6 * std::sin(-1.6);

                p["nonleptonic::Re{alpha4_c}@QCDF"] = 1.7 * std::cos(1.7);
                p["nonleptonic::Im{alpha4_c}@QCDF"] = 1.7 * std::sin(1.7);
                p["nonleptonic::Re{alpha3_c}@QCDF"] = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{alpha3_c}@QCDF"] = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{b4_c}@QCDF"]     = 1.9 * std::cos(1.9);
                p["nonleptonic::Im{b4_c}@QCDF"]     = 1.9 * std::sin(1.9);
                p["nonleptonic::Re{bS4_c}@QCDF"]    = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{bS4_c}@QCDF"]    = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                  = 0.0;


                static const double eps = 0.1e-6;

                Options o{
                    {     "q",    "u" },
                    {    "P1", "pi^0" },
                    {    "P2", "pi^+" },
                    { "model",  "CKM" },
                };

                QCDFRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(1.7700615877797436e-7, -3.729986432372129e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(1.824413161575804e-8, 2.9651504220674708e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(), complex<double>(1.9525029039373238e-7, -7.648360103046585e-9), eps);

                Options o2{
                    {     "q",    "d" },
                    {    "P1", "pi^+" },
                    {    "P2", "pi^-" },
                    { "model",  "CKM" },
                };

                QCDFRepresentation<PToPP> d2(p, o2);

                TEST_CHECK_RELATIVE_ERROR_C(d2.ordered_amplitude(), complex<double>(9.28644842045893e-10, 1.7567955154889877e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.inverse_amplitude(), complex<double>(2.5029133262201755e-7, -5.25160783494036e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.amplitude(), complex<double>(2.5121997746406343e-7, -5.23403987978547e-8), eps);


                Options o3{
                    {     "q",    "d" },
                    {    "P1", "pi^0" },
                    {    "P2", "pi^0" },
                    { "model",  "CKM" },
                };

                QCDFRepresentation<PToPP> d3(p, o3);

                TEST_CHECK_RELATIVE_ERROR_C(d3.ordered_amplitude(), complex<double>(-1.24690300570092e-8, -2.0761356952422998e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.inverse_amplitude(), complex<double>(-1.24690300570092e-8, -2.0761356952422998e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.amplitude(), complex<double>(-2.49380601140184e-8, -4.1522713904845995e-8), eps);

                Options o4{
                    {     "q",      "d" },
                    {    "P1",    "K_d" },
                    {    "P2", "Kbar_d" },
                    { "model",    "CKM" },
                };

                QCDFRepresentation<PToPP> d4(p, o4);

                TEST_CHECK_RELATIVE_ERROR_C(d4.ordered_amplitude(), complex<double>(1.6400147304326383e-7, 1.524228015749374e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.inverse_amplitude(), complex<double>(5.463352288201433e-10, 1.6014652984753204e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.amplitude(), complex<double>(1.6454780827208397e-7, 1.540242668734127e-8), eps);

                Options o5{
                    {            "q",      "s" },
                    {           "P1",    "K_d" },
                    {           "P2", "Kbar_d" },
                    {        "model",    "CKM" },
                    { "cp-conjugate",  "false" }
                };

                QCDFRepresentation<PToPP> d5(p, o5);

                TEST_CHECK_RELATIVE_ERROR_C(d5.ordered_amplitude(), complex<double>(-3.4371121163743705e-9, -1.1816622077486445e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.inverse_amplitude(), complex<double>(-7.822191406502231e-7, -1.036958299899132e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.amplitude(), complex<double>(-7.856562527665976e-7, -1.0487749219766184e-7), eps);


                Options o6{
                    {     "q",      "s" },
                    {    "P1", "Kbar_d" },
                    {    "P2",    "eta" },
                    { "model",    "CKM" },
                };

                QCDFRepresentation<PToPP> d6(p, o6);


                TEST_CHECK_RELATIVE_ERROR_C(d6.ordered_amplitude(), complex<double>(1.0387328442424082e-7, 2.4381522005734528e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.inverse_amplitude(), complex<double>(-8.803841392611456e-8, -8.123373428047469e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.amplitude(), complex<double>(1.5834870498126266e-8, 1.625814857768706e-8), eps);

                Options o7{
                    {     "q",         "d" },
                    {    "P1",       "K_d" },
                    {    "P2", "eta_prime" },
                    { "model",       "CKM" },
                };

                QCDFRepresentation<PToPP> d7(p, o7);


                TEST_CHECK_RELATIVE_ERROR_C(d7.ordered_amplitude(), complex<double>(-7.836971891724648e-7, 2.7757074325068556e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d7.inverse_amplitude(), complex<double>(-1.886061740105515e-7, -2.4826978870673098e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d7.amplitude(), complex<double>(-9.723033631830163e-7, 2.930095454395456e-9), eps);
            }
        }
} qcdf_amplitudes_test;
