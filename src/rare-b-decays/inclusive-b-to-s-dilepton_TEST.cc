/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/inclusive-b-to-s-dilepton.hh>

#include <cmath>
#include <limits>

using namespace test;
using namespace wf;

class BToXsDileptonLargeRecoilTest :
    public TestCase
{
    public:
        BToXsDileptonLargeRecoilTest() :
            TestCase("b_to_x_s_dilepton_HLMW2005_test")
        {
        }

        virtual void run() const
        {
            /* HLMW2005 */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["Re{c7}"] = -0.331;
                p["c8"] = -0.181;
                p["Re{c9}"] = +4.27;
                p["Re{c10}"] = -4.173;

                ObservableOptions oo;

                BToXsDilepton<HLMW2005> d(p, oo);

                const double eps = 1e-11;

                // Without log-enhanced em corrections: BR = 1.55709e-6
                TEST_CHECK_NEARLY_EQUAL(1.31892e-6, d.integrated_branching_ratio(1.00, 6.00), eps);
            }
        }
} b_to_x_s_dilepton_large_recoil_test;
