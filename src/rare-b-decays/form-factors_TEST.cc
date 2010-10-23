/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
