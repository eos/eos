/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2014 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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
#include <eos/form-factors/form-factors.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace eos;

class BCL2008FormFactorsTest :
    public TestCase
{
    public:
        BCL2008FormFactorsTest() :
            TestCase("bcl2008_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi@BCL2008", p);

                p["B->pi::f_+(0)@BCL2008"] = 1.0;
                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16483, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.40109, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.77364, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.47348, ff->f_p(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.17917, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.45663, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.94890, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.06978, ff->f_p(20.0), eps);
            }
        }
} bcl2008_form_factors_test;

class BToKstarBZ2004FormFactorsTest :
    public TestCase
{
    public:
        BToKstarBZ2004FormFactorsTest() :
            TestCase("b_to_kstar_bz2004_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*@BZ2004", p);
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.435248700565, ff->v( 1.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.460622528037, ff->v( 2.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.528633677061, ff->v( 4.3),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.589668921499, ff->v( 6.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.630935951233, ff->v( 7.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(1.113394334360, ff->v(14.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(1.367625759800, ff->v(16.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(2.034040927130, ff->v(19.2),   eps);

            TEST_CHECK_NEARLY_EQUAL(0.397077979401, ff->a_0( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.422486830218, ff->a_0( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.491687825118, ff->a_0( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.555085117357, ff->a_0( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.598608222322, ff->a_0( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.141522672280, ff->a_0(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.449013067940, ff->a_0(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(2.310470446100, ff->a_0(19.2), eps);

            TEST_CHECK_NEARLY_EQUAL(0.297364144236, ff->a_1( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.305112037520, ff->a_1( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.324562084257, ff->a_1( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.340610820244, ff->a_1( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.350814859197, ff->a_1( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.443904473086, ff->a_1(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.480319934372, ff->a_1(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.552889518414, ff->a_1(19.2), eps);

            TEST_CHECK_NEARLY_EQUAL(0.269896193772, ff->a_2( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.282547200000, ff->a_2( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.314867291642, ff->a_2( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.342079395085, ff->a_2( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.359608888889, ff->a_2( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.525473684211, ff->a_2(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.592222222222, ff->a_2(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.726406900654, ff->a_2(19.2), eps);
        }
} b_to_kstar_bz2004_form_factors_test;

class BsToPhiBZ2004FormFactorsTest :
    public TestCase
{
    public:
        BsToPhiBZ2004FormFactorsTest() :
            TestCase("bs_to_phi_bz2004_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("Bs->phi@BZ2004", p);
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.460064372742, ff->v( 1.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.487497702486, ff->v( 2.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.561398220530, ff->v( 4.3),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.628128476595, ff->v( 6.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(0.673439601066, ff->v( 7.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(1.210691300820, ff->v(14.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(1.496461092910, ff->v(16.0),   eps);
            TEST_CHECK_NEARLY_EQUAL(2.243708839070, ff->v(19.2),   eps);

            TEST_CHECK_NEARLY_EQUAL(0.501168940104, ff->a_0( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.529926892997, ff->a_0( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.608035423431, ff->a_0( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.679412484948, ff->a_0( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.728355734841, ff->a_0( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.339329806120, ff->a_0(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.687311851890, ff->a_0(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(2.669328444140, ff->a_0(19.2), eps);

            TEST_CHECK_NEARLY_EQUAL(0.316666291503, ff->a_1( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.325834394904, ff->a_1( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.349079404467, ff->a_1( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.368510805501, ff->a_1( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.380985781991, ff->a_1( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.499304347826, ff->a_1(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.547922103213, ff->a_1(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.649038062284, ff->a_1(19.2), eps);

            TEST_CHECK_NEARLY_EQUAL(0.245013923850, ff->a_2( 1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.256763995920, ff->a_2( 2.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.286954532316, ff->a_2( 4.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.312562021204, ff->a_2( 6.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.329147369735, ff->a_2( 7.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.489396953286, ff->a_2(14.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.555501255802, ff->a_2(16.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.691037087622, ff->a_2(19.2), eps);
        }
} bs_to_phi_bz2004_form_factors_test;


class BToKBZ2004FormFactorsTest :
    public TestCase
{
    public:
        BToKBZ2004FormFactorsTest() :
            TestCase("b_to_k_bz2004_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K@BZ2004v2", p);
            TEST_CHECK(0 != ff.get());

            // compare with results from mathematica
            TEST_CHECK_NEARLY_EQUAL(0.3795836109732980, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4357453886946407, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5081674405207227, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6042457948942926, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7362617433349268, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9259154765170850, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.2147208262068898, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.6909359437761546, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.5720556190766133, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.4079444647992955, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4697520247888276, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5496616371766049, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6559634649021799, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8024484632787752, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(1.0135299315425548, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.3360026902428106, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.8695529404263398, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.8603555237203935, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3515870307167236, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3761959829580037, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4045091623036650, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4374309978768577, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4761864406779661, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5224767540152155, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5787359550561798, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6485729275970619, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7375775656324581, ff->f_0(20.7), eps);
        }
} b_to_k_bz2004_form_factors_test;


class BToKBZ2004SplitFormFactorsTest :
    public TestCase
{
    public:
        BToKBZ2004SplitFormFactorsTest() :
            TestCase("b_to_k_bz2004_split_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K@BZ2004v2Split", p);
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.326792, ff->f_p( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.368753, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.421625, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.488304, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.574754, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.692499, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.863899, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.13432,  ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.60395,  ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.52611,  ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.349555, ff->f_t( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.395537, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.453490, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.526911, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.622869, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.754937, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.949382, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.25949,  ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.80311,  ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.87884,  ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.327521, ff->f_0( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.345076, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.364630, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.387123, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.414010, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.447569, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.491474, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.551925, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.640084, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.777994, ff->f_0(20.7), eps);
        }
} b_to_k_bz2004_split_form_factors_test;

class BToKstarKMPW2010FormFactorsTest:
    public TestCase
{
    public:
        BToKstarKMPW2010FormFactorsTest() :
            TestCase("b_to_kstar_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            // central values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.36; p["B->K^*::b^V_1@KMPW2010"]  =  -4.8;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.29; p["B->K^*::b^A0_1@KMPW2010"] = -18.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.25; p["B->K^*::b^A1_1@KMPW2010"] =  +0.34;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.23; p["B->K^*::b^A2_1@KMPW2010"] =  -0.85;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*@KMPW2010", p);
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.421964, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.497285, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.590298, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.707366, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.858208, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.058404, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.334545, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.735915, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.410875, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.560285, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.747894, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.988041, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.302880, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.728452, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.327350, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(3.218300, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.266631, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.285779, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.308075, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334383, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.365917, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.404440, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.452618, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514681, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.250198, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.273759, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.301563, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334822, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.375252, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.425358, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.488950, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572108, ff->a_2(18.4), eps);
            }

            // raised values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.59; p["B->K^*::b^V_1@KMPW2010"]  =  -4.0;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.39; p["B->K^*::b^A0_1@KMPW2010"] = -16.9;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.41; p["B->K^*::b^A1_1@KMPW2010"] =  +1.2;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.42; p["B->K^*::b^A2_1@KMPW2010"] =  +2.03;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*@KMPW2010", p);
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.683006, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.795816, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.934831, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.109443, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.333989, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.631445, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.041006, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.635294, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.543363, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.732820, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.970583, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.274770, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.673378, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(2.211932, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.969510, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(4.096084, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.430965, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.454709, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.481882, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.513361, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.550366, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.594652, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.648841, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.717050, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.435238, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.451990, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.470539, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.491252, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514616, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.541305, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572274, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.608954, ff->a_2(18.4), eps);
            }

            // lowered values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.24; p["B->K^*::b^V_1@KMPW2010"]  =  -5.2;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.22; p["B->K^*::b^A0_1@KMPW2010"] = -21.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.15; p["B->K^*::b^A1_1@KMPW2010"] =  -0.46;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.13; p["B->K^*::b^A2_1@KMPW2010"] =  -2.2;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*@KMPW2010", p);
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.2830474, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3354240, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4001625, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4817170, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.5868880, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.7265850, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.9194250, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.1999232, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.3236650, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4519470, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.6131980, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.8198090, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.0909360, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.4577380, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.9743500, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.7434500, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.1621258, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1762210, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1927954, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2125507, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2364774, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2660210, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3033756, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3520390, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.1445568, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1616854, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1820760, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2066815, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2368550, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2745790, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3228774, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3865910, ff->a_2(18.4), eps);
            }
        }
} b_to_kstar_kmpw2010_form_factors_test;

class BToKKMPW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKKMPW2010FormFactorsTest() :
            TestCase("b_to_k_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K@KMPW2010", p);
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.3844144474375, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4381960494587, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5043916755865, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5874824689246, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6943252577904, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8359400978914, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.0311821367700, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.3150955485260, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(1.7607545970800, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.4417893713152, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5045402114062, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5818230155226, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6788902047561, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8037805019019, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9694161202369, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.1979110863960, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.5303740143290, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.0525389704830, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3691335300041, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4001522415666, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4333067770539, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4689025169045, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5073172894003, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5490271105427, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5946449022903, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6449812193912, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7011445924499, ff->f_0(20.7), eps);
        }
} b_to_k_kmpw2010_form_factors_test;


class BToBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToBFW2010FormFactorsTest() :
            TestCase("b_to_k_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K@BFW2010", p);
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.3285329106071, ff->f_p( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3765270351772, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4352582092738, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5083556849771, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6012033334785, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7221134890713, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8845704673796, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.1118693561113, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.4477383401655, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(1.9842848427703, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3529456128368, ff->f_t( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4057232066867, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4703717236260, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5509131505349, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6533164641168, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7868011822009, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9663293753277, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.2177594231842, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.5896500926049, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.1843219169364, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3208496768629, ff->f_0( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3456358228669, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3724574045224, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4016438048517, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4336098586343, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4688881948925, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5081785681054, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5524262352940, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6029528628018, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6616893406433, ff->f_0(20.7), eps);
        }
} b_to_k_bfw2010_form_factors_test;

