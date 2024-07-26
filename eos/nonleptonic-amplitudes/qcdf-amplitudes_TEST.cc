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
#include <eos/nonleptonic-amplitudes/qcdf-amplitudes.hh>

using namespace test;
using namespace eos;

class QCDFAmplitudesTest :
    public TestCase
{
    public:
        QCDFAmplitudesTest() :
            TestCase("qcdf_amplitudes_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();

                p["nonleptonic::Re{alpha1}@QCDF"]     =  0.1 * std::cos( 0.1);
                p["nonleptonic::Im{alpha1}@QCDF"]     =  0.1 * std::sin( 0.1);
                p["nonleptonic::Re{alpha2}@QCDF"]     = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{alpha2}@QCDF"]     = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{b2}@QCDF"]         =  0.3 * std::cos( 0.3);
                p["nonleptonic::Im{b2}@QCDF"]         =  0.3 * std::sin( 0.3);
                p["nonleptonic::Re{b1}@QCDF"]         = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{b1}@QCDF"]         = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{bS2}@QCDF"]        =  0.5 * std::cos( 0.5);
                p["nonleptonic::Im{bS2}@QCDF"]        =  0.5 * std::sin( 0.5);
                p["nonleptonic::Re{bS1}@QCDF"]        = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{bS1}@QCDF"]        = -0.6 * std::sin(-0.6);

                p["nonleptonic::Re{alpha4_u}@QCDF"]   =  0.7 * std::cos( 0.7);
                p["nonleptonic::Im{alpha4_u}@QCDF"]   =  0.7 * std::sin( 0.7);
                p["nonleptonic::Re{alpha3_u}@QCDF"]   = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{alpha3_u}@QCDF"]   = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{b4_u}@QCDF"]       =  0.9 * std::cos( 0.9);
                p["nonleptonic::Im{b4_u}@QCDF"]       =  0.9 * std::sin( 0.9);  
                p["nonleptonic::Re{bS4_u}@QCDF"]      = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{bS4_u}@QCDF"]      = -1.0 * std::sin(-1.0); 

                p["nonleptonic::Re{alpha4EW_c}@QCDF"] =  1.1 * std::cos( 1.1);
                p["nonleptonic::Im{alpha4EW_c}@QCDF"] =  1.1 * std::sin( 1.1);    
                p["nonleptonic::Re{alpha3EW_c}@QCDF"] = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{alpha3EW_c}@QCDF"] = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{b3EW_c}@QCDF"]     =  1.3 * std::cos( 1.3);
                p["nonleptonic::Im{b3EW_c}@QCDF"]     =  1.3 * std::sin( 1.3);
                p["nonleptonic::Re{b4EW_c}@QCDF"]     = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{b4EW_c}@QCDF"]     = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{bS3EW_c}@QCDF"]    =  1.5 * std::cos( 1.5);
                p["nonleptonic::Im{bS3EW_c}@QCDF"]    =  1.5 * std::sin( 1.5);
                p["nonleptonic::Re{bS4EW_c}@QCDF"]    = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{bS4EW_c}@QCDF"]    = -1.6 * std::sin(-1.6); 
                
                p["nonleptonic::Re{alpha4_c}@QCDF"]   =  1.7 * std::cos( 1.7);
                p["nonleptonic::Im{alpha4_c}@QCDF"]   =  1.7 * std::sin( 1.7);                          
                p["nonleptonic::Re{alpha3_c}@QCDF"]   = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{alpha3_c}@QCDF"]   = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{b4_c}@QCDF"]       =  1.9 * std::cos( 1.9);
                p["nonleptonic::Im{b4_c}@QCDF"]       =  1.9 * std::sin( 1.9);
                p["nonleptonic::Re{bS4_c}@QCDF"]      = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{bS4_c}@QCDF"]      = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                    =  0.0;


                static const double eps = 0.1e-6;

                Options o
                {
                    { "q", "u" },
                    { "P1", "pi^0" },
                    { "P2", "pi^+" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d(p, o);

                TEST_CHECK_RELATIVE_ERROR_C(d.ordered_amplitude() , complex<double>( 1.7689401198311896e-7, -3.727624792246323e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.inverse_amplitude() , complex<double>( 1.8233115998322563e-8, 2.96339568973656e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d.amplitude() , complex<double>( 1.9512712798144152e-7, -7.642291025097634e-9), eps);

                Options o2
                {
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d2(p, o2);
                
                TEST_CHECK_RELATIVE_ERROR_C(d2.ordered_amplitude() , complex<double>(9.28644842045893e-10,1.7567955154889877e-10 ), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.inverse_amplitude() , complex<double>( 2.5013271447745307e-7, -5.248267579203517e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d2.amplitude() , complex<double>( 2.5106135931949895e-7, -5.230699624048627e-8), eps);


                Options o3
                {
                    { "q", "d" },
                    { "P1", "pi^0" },
                    { "P2", "pi^0" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d3(p, o3);
                
                TEST_CHECK_RELATIVE_ERROR_C(d3.ordered_amplitude() , complex<double>(-1.246123992443743e-8 ,-2.0748947663914343e-8 ), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.inverse_amplitude() , complex<double>( -1.246123992443743e-8, -2.0748947663914343e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d3.amplitude() , complex<double>( -2.492247984887486e-8, -4.1497895327828685e-8), eps);

                Options o4
                {
                    { "q", "d" },
                    { "P1", "K_d" },
                    { "P2", "Kbar_d" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d4(p, o4);
                
                TEST_CHECK_RELATIVE_ERROR_C(d4.ordered_amplitude() , complex<double>(5.463352288201433e-10 ,1.6014652984753204e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.inverse_amplitude() , complex<double>( 5.463352288201433e-10, 1.6014652984753204e-10), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d4.amplitude() , complex<double>( 1.0926704576402867e-9, 3.202930596950641e-10), eps);

                Options o5
                {
                    { "q", "s" },
                    { "P1", "K_d" },
                    { "P2", "Kbar_d" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d5(p, o5);

                TEST_CHECK_RELATIVE_ERROR_C(d5.ordered_amplitude() , complex<double>(-3.4371121163743705e-9 ,-1.1816622077486445e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.inverse_amplitude() , complex<double>( -7.761352825560144e-7, -1.0289498761479221e-7), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d5.amplitude() , complex<double>( -7.795723946723889e-7, -1.0407664982254085e-7), eps);


                Options o6
                {
                    { "q", "s" },
                    { "P1", "Kbar_d" },
                    { "P2", "eta" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                QCDFRepresentation<PToPP> d6(p, o6);

                
                TEST_CHECK_RELATIVE_ERROR_C(d6.ordered_amplitude() , complex<double>( 1.0288965971149677e-7, 2.415064197039064e-8), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.inverse_amplitude() , complex<double>(-8.735065637966213e-8 ,-8.059913500401798e-9), eps);
                TEST_CHECK_RELATIVE_ERROR_C(d6.amplitude() , complex<double>( 1.553900333183462e-8, 1.6090728469988845e-8), eps);

            }
        }
} qcdf_amplitudes_test;

