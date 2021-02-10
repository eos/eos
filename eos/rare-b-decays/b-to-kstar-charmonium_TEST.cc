#include <test/test.hh>
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BToKstarCharmoniumGvDV2020Test :
    public TestCase
{
    public:
    BToKstarCharmoniumGvDV2020Test() :
            TestCase("b_to_kstar_charmonium_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]                               = 5.27942;
            p["mass::K_d^*"]                             = 0.89555;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 9.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

            p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 2.0;
            p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 3.0;
            p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 4.0;
            p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 5.0;
            p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 6.0;
            p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 7.0;
            p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 8.0;
            p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 9.0;
            p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 10.0;
            p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 11.0;
            p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 12.0;
            p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 13.0;
            p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 14.0;
            p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 15.0;
            p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 16.0;
            p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 17.0;
            p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 18.0;
            p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 19.0;

            Options oo;
            oo.set("model",               "WilsonScan");
            oo.set("q",                   "d");
            oo.set("nonlocal-formfactor", "GvDV2020");
            oo.set("psi",                 "J/psi");

            BToKstarCharmonium c(p, oo);

            TEST_CHECK_NEARLY_EQUAL(c.branching_ratio(),  12481750, 1.);

            TEST_CHECK_NEARLY_EQUAL(c.S_1c_LHCb(),  0.18679,    eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_1s_LHCb(),  0.60990,    eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_3_LHCb(),  -0.23336,    eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_4_LHCb(),   0.24446,    eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_8_LHCb(),  -0.00635388, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_9_LHCb(),   0.00870404, eps);

        }
} b_to_kstar_charmonium_GvDV2020_test;


class BToKstarCharmoniumGRvDV2021Test :
    public TestCase
{
    public:
    BToKstarCharmoniumGRvDV2021Test() :
            TestCase("b_to_kstar_charmonium_GRvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]                                = 5.27942;
            p["mass::K_d^*"]                              = 0.89555;
            p["mass::J/psi"]                              = 3.0969;
            p["mass::psi(2S)"]                            = 3.6860;
            p["mass::D^0"]                                = 1.86723;
            p["b->sccbar::t_0"]                           = 9.0;
            p["b->sccbar::t_s"]                           = -17.4724;
            p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;

            p["B->K^*ccbar::Re{alpha_0^perp}@GRvDV2021"]  = 2.0;
            p["B->K^*ccbar::Im{alpha_0^perp}@GRvDV2021"]  = 3.0;
            p["B->K^*ccbar::Re{alpha_1^perp}@GRvDV2021"]  = 4.0;
            p["B->K^*ccbar::Im{alpha_1^perp}@GRvDV2021"]  = 5.0;
            p["B->K^*ccbar::Re{alpha_2^perp}@GRvDV2021"]  = 6.0;
            p["B->K^*ccbar::Im{alpha_2^perp}@GRvDV2021"]  = 7.0;
            p["B->K^*ccbar::Re{alpha_0^para}@GRvDV2021"]  = 8.0;
            p["B->K^*ccbar::Im{alpha_0^para}@GRvDV2021"]  = 9.0;
            p["B->K^*ccbar::Re{alpha_1^para}@GRvDV2021"]  = 10.0;
            p["B->K^*ccbar::Im{alpha_1^para}@GRvDV2021"]  = 11.0;
            p["B->K^*ccbar::Re{alpha_2^para}@GRvDV2021"]  = 12.0;
            p["B->K^*ccbar::Im{alpha_2^para}@GRvDV2021"]  = 13.0;
            p["B->K^*ccbar::Re{alpha_0^long}@GRvDV2021"]  = 14.0;
            p["B->K^*ccbar::Im{alpha_0^long}@GRvDV2021"]  = 15.0;
            p["B->K^*ccbar::Re{alpha_1^long}@GRvDV2021"]  = 16.0;
            p["B->K^*ccbar::Im{alpha_1^long}@GRvDV2021"]  = 17.0;
            p["B->K^*ccbar::Re{alpha_2^long}@GRvDV2021"]  = 18.0;
            p["B->K^*ccbar::Im{alpha_2^long}@GRvDV2021"]  = 19.0;

            Options oo;
            oo.set("model",               "WilsonScan");
            oo.set("q",                   "d");
            oo.set("nonlocal-formfactor", "GRvDV2021");
            oo.set("psi",                 "J/psi");

            BToKstarCharmonium c(p, oo);

            TEST_CHECK_NEARLY_EQUAL(c.branching_ratio(),  890412., 1.);

            TEST_CHECK_NEARLY_EQUAL(c.S_1c_LHCb(),  0.263661, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_1s_LHCb(),  0.552254, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_3_LHCb(),  -0.309413, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_4_LHCb(),   0.298785, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_8_LHCb(),  -0.014683, eps);
            TEST_CHECK_NEARLY_EQUAL(c.S_9_LHCb(),   0.020114, eps);

        }
} b_to_kstar_charmonium_GRvDV2021_test;