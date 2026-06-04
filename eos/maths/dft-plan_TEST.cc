/*
 * Copyright (c) 2026 Danny van Dyk
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
#include <eos/maths/dft-plan-impl.hh>


using namespace test;
using namespace eos;

class DftPlanImplementationDetailsTestCase :
    public TestCase
{
    public:
        DftPlanImplementationDetailsTestCase() :
            TestCase("dft_plan_implementation_details_test")
        {
        }

        virtual void run() const override
        {
            // 1D
            {
                std::array<std::size_t, 1> time_domain_dimensions = { 16 };

                std::array<std::size_t, 1> frequency_domain_dimensions = dft::impl::frequency_dimensions<1>(time_domain_dimensions);

                TEST_CHECK(frequency_domain_dimensions[0] == 9);
            }

            // 2D
            {
                std::array<std::size_t, 2> time_domain_dimensions = { 32, 64 };

                std::array<std::size_t, 2> frequency_domain_dimensions = dft::impl::frequency_dimensions<2>(time_domain_dimensions);

                TEST_CHECK(frequency_domain_dimensions[0] == 32);
                TEST_CHECK(frequency_domain_dimensions[1] == 33);
            }
        }
} dft_plan_implementation_details_test;

class DftPlanOddDimensionsTest :
    public TestCase
{
    public:
        DftPlanOddDimensionsTest() :
            TestCase("dft_plan_odd_dimensions_test")
        {
        }

        virtual void run() const override
        {
            // odd dimensions must throw upon plan construction

            // 1D
            {
                std::array<std::size_t, 1> dimensions{ 15 };
                TEST_CHECK_THROWS(eos::InternalError, (dft::Plan<1, dft::Direction::Forward>(dimensions)));
            }

            // 2D: odd in the last dimension
            {
                std::array<std::size_t, 2> dimensions{ 4, 7 };
                TEST_CHECK_THROWS(eos::InternalError, (dft::Plan<2, dft::Direction::Forward>(dimensions)));
            }

            // 2D: odd in a leading dimension
            {
                std::array<std::size_t, 2> dimensions{ 3, 8 };
                TEST_CHECK_THROWS(eos::InternalError, (dft::Plan<2, dft::Direction::Backward>(dimensions)));
            }

            // 4D: odd in an interior dimension
            {
                std::array<std::size_t, 4> dimensions{ 2, 4, 5, 4 };
                TEST_CHECK_THROWS(eos::InternalError, (dft::Plan<4, dft::Direction::Forward>(dimensions)));
            }
        }
} dft_plan_odd_dimensions_test;

class DftPlanRankOneTest :
    public TestCase
{
    public:
        DftPlanRankOneTest() :
            TestCase("dft_plan_rank_one_test")
        {
        }

        virtual void run() const override
        {
            static const std::array<std::size_t, 1> dimensions{ 256 };
            static const std::array<double, 256> data{
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.00390244,
                 0.0310059,  0.0577354,  0.0840912,  0.110073,   0.135681,   0.160915,   0.185776,   0.210262,
                 0.234375,   0.258114,   0.281479,   0.30447,    0.327087,   0.349331,   0.371201,   0.392696,
                 0.413818,   0.434566,   0.454941,   0.474941,   0.494568,   0.513821,   0.5327,     0.551205,
                 0.569336,   0.587093,   0.604477,   0.621487,   0.638123,   0.654385,   0.670273,   0.685787,
                 0.700928,   0.715694,   0.730087,   0.744106,   0.75856,    0.771854,   0.784775,   0.797323,
                 0.809499,   0.821302,   0.832732,   0.84379,    0.854476,   0.864789,   0.874729,   0.884298,
                 0.893494,   0.902319,   0.910772,   0.918853,   0.926562,   0.9339,     0.940867,   0.947463,
                 0.953689,   0.959544,   0.965028,   0.970143,   0.974888,   0.979264,   0.983271,   0.986909,
                 0.99018,    0.993083,   0.99562,    0.99779,    0.999595,   1.00104,    1.00211,    1.00283,
                 1.00318,    1.00318,    1.00281,    1.00209,    1.00102,    0.999598,   0.997825,   0.995709,
                 0.993251,   0.990459,   0.987337,   0.983892,   0.980133,   0.976071,   0.971718,   0.967089,
                 0.962205,   0.957091,   0.951777,   0.946304,   0.940727,   0.935117,   0.929569,   0.924217,
                 0.919249,   0.914938,   0.911691,   0.910138,   0.911296,   0.916891,   0.930053,   0.956933,
                 1.01102,    1.12675,    1.41501,    2.40647,    10.0072,    14.1563,    2.68786,    1.42151,
                 1.05902,    0.90139,    0.814085,   0.757068,   0.715084,   0.681302,   0.652295,   0.626191,
                 0.601898,   0.578742,   0.556291,   0.534257,   0.512443,   0.490707,   0.46895,    0.447097,
                 0.425093,   0.402894,   0.380469,   0.357792,   0.334842,   0.311603,   0.288062,   0.264208,
                 0.240033,   0.215529,   0.190691,   0.165513,   0.139991,   0.114121,   0.0879005,  0.0613267,
                 0.0343973,  0.00711023, 0.00303869, 0.00288261, 0.00273826, 0.00260448, 0.00248028, 0.00236475,
                 0.00225711, 0.00215665, 0.00206276, 0.00197486, 0.00189247, 0.00181512, 0.00174242, 0.00167401,
                 0.00160954, 0.00154873, 0.00149131, 0.00143701, 0.00138563, 0.00133696, 0.0012908,  0.001247,
                 0.00120539, 0.00116583, 0.00112818, 0.00109233, 0.00105816, 0.00102556, 0.000994456,0.000964742,
                 0.000936339,0.000909173,0.000883172,0.000858271,0.000834408,0.000811527,0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,
                 0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0
            };

            // Create plan for forward transformation
            dft::Plan<1, dft::Direction::Forward> forward_plan(dimensions);
            dft::Container<double, 1> & time_domain_container = forward_plan.time_domain_container();
            dft::Container<std::complex<double>, 1> & frequency_domain_container = forward_plan.frequency_domain_container();

            TEST_CHECK_EQUAL(dimensions, time_domain_container.dimensions());
            TEST_CHECK_EQUAL(dimensions[0] / 2 + 1, frequency_domain_container.dimensions()[0]);

            // Execute forward plan
            std::copy(data.begin(), data.end(), forward_plan.time_domain_container().data());
            forward_plan.transform();

            // Expected values obtained from Mathematica's Fourier[] with the FFTW3 convention
            // (negative exponent, no normalisation).
            static const std::array<std::complex<double>, 129> expected{{
                {  127.799,       0.0        },  // k=0
                {  -88.256,      23.5526     },  // k=1
                {   11.0419,    -28.9844     },  // k=2
                {   34.2744,     12.6913     },  // k=3
                {  -29.1046,     12.7161     },  // k=4
                {    2.84081,   -28.0943     },  // k=5
                {   18.1886,     22.1584     },  // k=6
                {  -26.2979,      0.16831    },  // k=7
                {   17.6264,    -21.8945     },  // k=8
                {    5.35677,    26.7939     },  // k=9
                {  -25.2374,    -11.5946     },  // k=10
                {   24.5858,    -11.9222     },  // k=11
                {   -5.06143,    26.0222     },  // k=12
                {  -16.7516,    -20.3811     },  // k=13
                {   25.5603,     -0.337285   },  // k=14
                {  -15.8668,     20.3451     },  // k=15
                {   -5.72062,   -24.7116     },  // k=16
                {   23.0344,     10.5391     },  // k=17
                {  -22.3772,     11.177      },  // k=18
                {    4.78953,   -24.0565     },  // k=19
                {   15.5613,     18.6817     },  // k=20
                {  -23.7924,      0.507591   },  // k=21
                {   14.4204,    -18.8779     },  // k=22
                {    5.60399,    22.7266     },  // k=23
                {  -21.2193,     -9.52464    },  // k=24
                {   20.4473,    -10.4762     },  // k=25
                {   -4.29418,    22.1864     },  // k=26
                {  -14.4654,    -17.0555     },  // k=27
                {   21.9519,     -0.679903   },  // k=28
                {  -13.0689,     17.4856     },  // k=29
                {   -5.41148,   -20.8311     },  // k=30
                {   19.549,       8.54759    },  // k=31
                {  -18.644,       9.8169     },  // k=32
                {    3.75532,   -20.4041     },  // k=33
                {   13.4365,     15.4963     },  // k=34
                {  -20.1502,      0.854899   },  // k=35
                {   11.7799,    -16.1627     },  // k=36
                {    5.20868,    19.0179     },  // k=37
                {  -17.9761,     -7.60421    },  // k=38
                {   16.9306,     -9.19644    },  // k=39
                {   -3.20933,    18.7025     },  // k=40
                {  -12.4645,    -13.9983     },  // k=41
                {   18.4092,     -1.03327    },  // k=42
                {  -10.5418,     14.9037     },  // k=43
                {   -5.01323,   -17.2796     },  // k=44
                {   16.4828,      6.69083    },  // k=45
                {  -15.2914,      8.61234    },  // k=46
                {    2.6662,    -17.0748     },  // k=47
                {   11.5436,     12.5555     },  // k=48
                {  -16.7321,      1.21572    },  // k=49
                {    9.34737,   -13.7035     },  // k=50
                {    4.83112,    15.6097     },  // k=51
                {  -15.0587,     -5.80385    },  // k=52
                {   13.7163,     -8.06227    },  // k=53
                {   -2.12857,    15.5146     },  // k=54
                {  -10.6693,    -11.1622     },  // k=55
                {   15.1165,     -1.40298    },  // k=56
                {   -8.1909,     12.5576     },  // k=57
                {   -4.66458,   -14.0013     },  // k=58
                {   13.6964,      4.93977    },  // k=59
                {  -12.1973,      7.54405    },  // k=60
                {    1.59653,   -14.0156     },  // k=61
                {    9.83774,     9.81292    },  // k=62
                {  -13.5579,      1.59577    },  // k=63
                {    7.06739,   -11.4612     },  // k=64
                {    4.51438,    12.4482     },  // k=65
                {  -12.3896,     -4.09519    },  // k=66
                {   10.7278,     -7.05563    },  // k=67
                {   -1.06908,    12.572      },  // k=68
                {   -9.04534,    -8.50242    },  // k=69
                {   12.0511,     -1.79487    },  // k=70
                {   -5.97217,    10.4102     },  // k=71
                {   -4.38071,   -10.9443     },  // k=72
                {   11.1325,      3.26677    },  // k=73
                {   -9.30152,     6.59508    },  // k=74
                {    0.544775,  -11.1781     },  // k=75
                {    8.2889,      7.22548    },  // k=76
                {  -10.5909,      2.00106    },  // k=77
                {    4.90077,    -9.40025    },  // k=78
                {    4.26349,     9.48357    },  // k=79
                {   -9.91992,    -2.45125    },  // k=80
                {    7.91245,    -6.16058    },  // k=81
                {   -0.021896,    9.8283     },  // k=82
                {   -7.56535,    -5.97707    },  // k=83
                {    9.17184,    -2.21514    },  // k=84
                {   -3.84889,     8.42744    },  // k=85
                {   -4.16252,    -8.0603     },  // k=86
                {    8.7469,      1.64541    },  // k=87
                {   -6.55501,     5.75041    },  // k=88
                {   -0.501393,   -8.51732    },  // k=89
                {    6.8718,      4.75226    },  // k=90
                {   -7.78847,     2.43798    },  // k=91
                {    2.8123,     -7.48792    },  // k=92
                {    4.07757,     6.66887    },  // k=93
                {   -7.60869,    -0.846065   },  // k=94
                {    5.22373,    -5.36296    },  // k=95
                {    1.02702,     7.23999    },  // k=96
                {   -6.20548,    -3.54621    },  // k=97
                {    6.43548,    -2.67044    },  // k=98
                {   -1.78689,     6.57797    },  // k=99
                {   -4.00841,    -5.30376    },  // k=100
                {    6.5007,      0.0500586  },  // k=101
                {   -3.91329,     4.99668    },  // k=102
                {   -1.55695,    -5.99125    },  // k=103
                {    5.56374,     2.35416    },  // k=104
                {   -5.10763,     2.91345    },  // k=105
                {    0.768571,   -5.694      },  // k=106
                {    3.95484,     3.9596     },  // k=107
                {   -5.41852,     0.745752   },  // k=108
                {    2.61847,    -4.65014    },  // k=109
                {    2.09323,     4.76618    },  // k=110
                {   -4.94405,    -1.1714     },  // k=111
                {    3.79975,    -3.16797    },  // k=112
                {    0.246689,    4.83253    },  // k=113
                {   -3.91671,    -2.63108    },  // k=114
                {    4.35784,    -1.54451    },  // k=115
                {   -1.33414,     4.32196    },  // k=116
                {   -2.63792,    -3.55992    },  // k=117
                {    4.34395,    -0.00673168 },  // k=118
                {   -2.5067,      3.435      },  // k=119
                {   -1.26292,    -3.99013    },  // k=120
                {    3.89387,     1.31295    },  // k=121
                {   -3.31442,     2.34937    },  // k=122
                {    0.0551965,  -4.01086    },  // k=123
                {    3.19315,     2.36773    },  // k=124
                {   -3.76109,     1.18489    },  // k=125
                {    1.2234,     -3.71559    },  // k=126
                {    2.28415,     3.1635     },  // k=127
                {   -3.88627,     0.0        },  // k=128
            }};

            for (std::size_t k = 0; k < 129; ++k)
            {
                TEST_CHECK_NEARLY_EQUAL(frequency_domain_container.data()[k], expected[k], 1.0e-2);
            }
        }
} dft_plan_rank_one_test;

class DftPlanRankTwoTest :
    public TestCase
{
    public:
        DftPlanRankTwoTest() :
            TestCase("dft_plan_rank_two_test")
        {
        }

        virtual void run() const override
        {
            static const std::array<std::size_t, 2> dimensions{ 4, 8 };

            // check dimension reporting
            {
                dft::Plan<2, dft::Direction::Forward> forward_plan(dimensions);

                TEST_CHECK_EQUAL(dimensions, forward_plan.time_domain_container().dimensions());

                static const std::array<std::size_t, 2> expected_freq_dims{ 4, 5 };
                TEST_CHECK_EQUAL(expected_freq_dims, forward_plan.frequency_domain_container().dimensions());
            }

            // round-trip test: forward followed by backward (with 1/N normalisation) recovers the input
            {
                static const std::array<double, 4 * 8> input{{
                     0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
                     8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
                    16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
                    24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0
                }};

                dft::Plan<2, dft::Direction::Forward> forward_plan(dimensions);
                std::copy(input.begin(), input.end(), forward_plan.time_domain_container().data());
                forward_plan.transform();

                dft::Plan<2, dft::Direction::Backward> backward_plan(dimensions);
                std::copy(
                    forward_plan.frequency_domain_container().data(),
                    forward_plan.frequency_domain_container().data() + 4 * 5,
                    backward_plan.frequency_domain_container().data());
                backward_plan.transform();

                const double norm = 1.0 / (4.0 * 8.0);
                const double eps = 1.0e-10;
                for (std::size_t n = 0; n < 4 * 8; ++n)
                {
                    TEST_CHECK_NEARLY_EQUAL(backward_plan.time_domain_container().data()[n] * norm, input[n], eps);
                }
            }
        }
} dft_plan_rank_two_test;

class DftPlanRankThreeTest :
    public TestCase
{
    public:
        DftPlanRankThreeTest() :
            TestCase("dft_plan_rank_three_test")
        {
        }

        virtual void run() const override
        {
            static const std::array<std::size_t, 3> dimensions{ 4, 4, 8 };

            // check dimension reporting
            {
                dft::Plan<3, dft::Direction::Forward> forward_plan(dimensions);

                TEST_CHECK_EQUAL(dimensions, forward_plan.time_domain_container().dimensions());

                static const std::array<std::size_t, 3> expected_freq_dims{ 4, 4, 5 };
                TEST_CHECK_EQUAL(expected_freq_dims, forward_plan.frequency_domain_container().dimensions());
            }

            // round-trip test: forward followed by backward (with 1/N normalisation) recovers the input
            {
                dft::Plan<3, dft::Direction::Forward> forward_plan(dimensions);
                double * time_in = forward_plan.time_domain_container().data();
                for (std::size_t n = 0; n < 4 * 4 * 8; ++n)
                {
                    time_in[n] = static_cast<double>(n);
                }

                forward_plan.transform();

                dft::Plan<3, dft::Direction::Backward> backward_plan(dimensions);
                std::copy(
                    forward_plan.frequency_domain_container().data(),
                    forward_plan.frequency_domain_container().data() + 4 * 4 * 5,
                    backward_plan.frequency_domain_container().data());
                backward_plan.transform();

                const double norm = 1.0 / (4.0 * 4.0 * 8.0);
                const double eps = 1.0e-10;
                for (std::size_t n = 0; n < 4 * 4 * 8; ++n)
                {
                    TEST_CHECK_NEARLY_EQUAL(backward_plan.time_domain_container().data()[n] * norm, static_cast<double>(n), eps);
                }
            }
        }
} dft_plan_rank_three_test;

class DftPlanRankFourTest :
    public TestCase
{
    public:
        DftPlanRankFourTest() :
            TestCase("dft_plan_rank_four_test")
        {
        }

        virtual void run() const override
        {
            static const std::array<std::size_t, 4> dimensions{ 2, 4, 4, 4 };

            // check dimension reporting
            {
                dft::Plan<4, dft::Direction::Forward> forward_plan(dimensions);

                TEST_CHECK_EQUAL(dimensions, forward_plan.time_domain_container().dimensions());

                static const std::array<std::size_t, 4> expected_freq_dims{ 2, 4, 4, 3 };
                TEST_CHECK_EQUAL(expected_freq_dims, forward_plan.frequency_domain_container().dimensions());
            }

            // round-trip test: forward followed by backward (with 1/N normalisation) recovers the input
            {
                dft::Plan<4, dft::Direction::Forward> forward_plan(dimensions);
                double * time_in = forward_plan.time_domain_container().data();
                for (std::size_t n = 0; n < 2 * 4 * 4 * 4; ++n)
                {
                    time_in[n] = static_cast<double>(n);
                }

                forward_plan.transform();

                dft::Plan<4, dft::Direction::Backward> backward_plan(dimensions);
                std::copy(
                    forward_plan.frequency_domain_container().data(),
                    forward_plan.frequency_domain_container().data() + 2 * 4 * 4 * 3,
                    backward_plan.frequency_domain_container().data());
                backward_plan.transform();

                const double norm = 1.0 / (2.0 * 4.0 * 4.0 * 4.0);
                const double eps = 1.0e-10;
                for (std::size_t n = 0; n < 2 * 4 * 4 * 4; ++n)
                {
                    TEST_CHECK_NEARLY_EQUAL(backward_plan.time_domain_container().data()[n] * norm, static_cast<double>(n), eps);
                }
            }
        }
} dft_plan_rank_four_test;
