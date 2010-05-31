/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace wf;

class BToKstarDileptonLargeRecoilTest :
    public TestCase
{
    public:
        BToKstarDileptonLargeRecoilTest() :
            TestCase("b_to_kstar_dilepton_large_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Large Recoil */
            {
                Parameters p = Parameters::Defaults();
                p["c7"] = -0.331;
                p["c8"] = -0.181;
                p["c9"] = +4.27;
                p["c10"] = -4.173;

                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LargeRecoil> d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(2.00, 4.30) * 1e7,      +1.1107, eps); // 1.1125 for m_b_PS = 4.6
                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), +0.0846, eps); // 0.0911 for m_b_PS = 4.6
                TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(2.00, 4.30),  +0.7833, eps); // 0.7800 for m_b_PS = 4.6
            }
        }
} b_to_kstar_dilepton_large_recoil_test;

class BToKstarDileptonLowRecoilTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil */
            {
                Parameters p = Parameters::Defaults();
                p["c7"] = -0.336;
                p["c9"] = +4.27;
                p["c10"] = -4.173;

                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), -0.4091, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(14.00, 19.21),  +0.3497, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(14.00, 19.21),     -0.4840, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(14.00, 19.21),     +1.6903, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(14.00, 19.21),     -0.5755, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(14.00, 19.21),                        +0.9967, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(14.00, 19.21),                        -0.9728, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(14.00, 19.21),                        -0.9588, eps);
            }
        }
} b_to_kstar_dilepton_low_recoil_test;
