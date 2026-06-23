/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
 * Copyright (c) 2025 Danny van Dyk
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
                /*p["eta::theta_18"]                  = 0.0;*/
                p["eta::theta_FKS"]                 = 0.686;


                static const double eps = 0.1e-6;

                Options o{
                    {     "q"_ok,    "u" },
                    {    "P1"_ok, "pi^0" },
                    {    "P2"_ok, "pi^+" },
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude(), complex<double>(1.7700615877797436e-7, -3.729986432372129e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude(), complex<double>(1.824413161575804e-8, 2.9651504220674708e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude(), complex<double>(1.9525029039373238e-7, -7.648360103046585e-9), eps);

                Options o2{
                    {     "q"_ok,    "d" },
                    {    "P1"_ok, "pi^+" },
                    {    "P2"_ok, "pi^-" },
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d2(p, o2);

                TEST_CHECK_RELATIVE_ERROR_C(d2.ordered_amplitude(), complex<double>(9.28644842045893e-10, 1.7567955154889877e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.inverse_amplitude(), complex<double>(2.5029133262201755e-7, -5.25160783494036e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.amplitude(), complex<double>(2.5121997746406343e-7, -5.23403987978547e-8), eps);


                Options o3{
                    {     "q"_ok,    "d" },
                    {    "P1"_ok, "pi^0" },
                    {    "P2"_ok, "pi^0" },
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d3(p, o3);

                TEST_CHECK_RELATIVE_ERROR_C(d3.ordered_amplitude(), complex<double>(-1.24690300570092e-8, -2.0761356952422998e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.inverse_amplitude(), complex<double>(-1.24690300570092e-8, -2.0761356952422998e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.amplitude(), complex<double>(-2.49380601140184e-8, -4.1522713904845995e-8), eps);

                Options o4{
                    {     "q"_ok,      "d" },
                    {    "P1"_ok,    "K_d" },
                    {    "P2"_ok, "Kbar_d" },
                    { "model"_ok,    "CKM" },
                };

                QCDFRepresentation<PToPP> d4(p, o4);

                TEST_CHECK_RELATIVE_ERROR_C(d4.ordered_amplitude(), complex<double>(1.6395931854535794e-7, 1.5238390528945357e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.inverse_amplitude(), complex<double>(5.463352288201433e-10, 1.60146529847532e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.amplitude(), complex<double>(1.645056537741781e-7, 1.539853705879289e-8), eps);

                Options o5{
                    {            "q"_ok,      "s" },
                    {           "P1"_ok,    "K_d" },
                    {           "P2"_ok, "Kbar_d" },
                    {        "model"_ok,    "CKM" },
                    { "cp-conjugate"_ok,  "false" }
                };

                QCDFRepresentation<PToPP> d5(p, o5);

                TEST_CHECK_RELATIVE_ERROR_C(d5.ordered_amplitude(), complex<double>(-3.4371121163743705e-9, -1.1816622077486445e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.inverse_amplitude(), complex<double>(-7.824200375723936e-7, -1.0372227484872056e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.amplitude(), complex<double>(-7.85857149688768e-7, -1.049039370564692e-7), eps);

                Options o6{
                    {     "q"_ok,    "u" },
                    {    "P1"_ok,  "eta" },
                    {    "P2"_ok, "pi^+" },
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d6(p, o6);

                TEST_CHECK_RELATIVE_ERROR_C(d6.ordered_amplitude(), complex<double>(1.046345494732723e-7, -2.205519360498948e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.inverse_amplitude(), complex<double>(2.1705544131419826e-7,1.8836775006222694e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.amplitude(), complex<double>(3.2168999078747063e-7,-3.2184185987667867e-9), eps);


                Options o7{
                    {     "q"_ok,    "s" },
                    {    "P1"_ok,  "Kbar_d" },
                    {    "P2"_ok,  "eta_prime"},
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d7(p, o7);


                TEST_CHECK_RELATIVE_ERROR_C(d7.ordered_amplitude(), complex<double>(4.753880811553182e-7,-9.062079140448238e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d7.inverse_amplitude(), complex<double>(2.32823512972189e-8,2.148280795119971e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d7.amplitude(), complex<double>(4.986704324525372e-7,-8.847251060936242e-8), eps);


                Options o8{
                    {     "q"_ok,    "d" },
                    {    "P1"_ok,  "eta" },
                    {    "P2"_ok, "eta_prime" },
                    { "model"_ok,  "CKM" },
                };

                QCDFRepresentation<PToPP> d8(p, o8);

                TEST_CHECK_RELATIVE_ERROR_C(d8.ordered_amplitude(), complex<double>(1.9773659594185027e-7,-3.7772480742497153e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d8.inverse_amplitude(), complex<double>(7.452553993234181e-8,6.2510416243630515e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d8.amplitude(), complex<double>(2.7226213587419205e-7,-3.15214391181341e-8), eps);

                Options o9
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta_prime" },
                    { "P2"_ok, "K_S" },
                };

                QCDFRepresentation<PToPP> d9(p, o9);

                TEST_CHECK_RELATIVE_ERROR_C(d9.ordered_amplitude(), complex<double>(-1.8653040020760233e-7,-2.455373653056162e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d9.inverse_amplitude(), complex<double>(-1.786584610102707e-6,7.126348192849077e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d9.amplitude(), complex<double>(-1.9731150103103095e-6,4.670974539792915e-8), eps);


            }
        }
} qcdf_amplitudes_test;
