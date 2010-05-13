/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cmath>
#include <limits>

using namespace test;
using namespace wf;

class BToKstarDileptonTest :
    public TestCase
{
    public:
        BToKstarDileptonTest() :
            TestCase("b_to_kstar_dilepton_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil */
            {
                Parameters p = Parameters::Defaults();
                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");
                BToKstarDilepton<LowRecoil> d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), -0.4087, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(14.00, 19.21),  +0.3497, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(14.00, 19.21),     -0.4839, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(14.00, 19.21),     +1.6902, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(14.00, 19.21),     -0.5750, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(14.00, 19.21),                        +0.9967, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(14.00, 19.21),                        -0.9718, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(14.00, 19.21),                        -0.9578, eps);
            }
        }
} b_to_kstar_dilepton_test;
