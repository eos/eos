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
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;

                ObservableOptions oo;

                BToXsDilepton<HLMW2005> d(p, oo);

                const double eps = 1e-11;

                TEST_CHECK_NEARLY_EQUAL(1.61522e-06, d.integrated_branching_ratio(1.00, 6.00), eps);
            }
        }
} b_to_x_s_dilepton_large_recoil_test;
