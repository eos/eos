/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <src/rare-b-decays/form-factors.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace eos;

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

