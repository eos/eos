/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/matrix.hh>

#include <cmath>
#include <array>

using namespace test;
using namespace eos;

class MatrixMultiplicationTest :
    public TestCase
{
    public:
        MatrixMultiplicationTest() :
            TestCase("matrix_multiplication_test")
        {
        }

        virtual void run() const
        {
            using std::array;

            // matrix times matrix
            {
                static const array<array<double, 4>, 6> x
                {{
                    array<double, 4>{{0.256065,   0.277201,   0.406745,   0.188430}},
                    array<double, 4>{{0.099332,   0.576983,   0.077084,   0.279125}},
                    array<double, 4>{{0.433378,   0.987969,   0.198432,   0.497537}},
                    array<double, 4>{{0.795675,   0.253245,   0.740202,   0.958491}},
                    array<double, 4>{{0.732521,   0.606855,   0.893422,   0.532790}},
                    array<double, 4>{{0.935207,   0.889252,   0.079741,   0.098048}},
                }};
                static const array<array<double, 3>, 4> y
                {{
                    array<double, 3>{{0.280209,   0.909076,   0.191122}},
                    array<double, 3>{{0.650764,   0.967410,   0.044342}},
                    array<double, 3>{{0.474950,   0.838596,   0.210205}},
                    array<double, 3>{{0.421354,   0.060448,   0.267996}},
                }};
                static const array<array<double, 3>, 6> z
                {{
                    array<double, 3>{{0.52472,   0.85343,   0.19723}},
                    array<double, 3>{{0.55753,   0.72999,   0.13558}},
                    array<double, 3>{{1.06826,   1.54622,   0.30169}},
                    array<double, 3>{{1.14318,   1.64699,   0.57577}},
                    array<double, 3>{{1.24900,   2.03442,   0.49750}},
                    array<double, 3>{{0.91993,   1.78324,   0.26121}},
                }};

                array<array<double, 3>, 6> result = x * y;

                for (unsigned i(0) ; i < 6 ; ++i)
                {
                    for (unsigned j(0) ; j < 3 ; ++j)
                    {
                        TEST_CHECK_NEARLY_EQUAL(result[i][j], z[i][j], 1e-5);
                    }
                }
            }

            // matrix times scalar
            {
                static const double x = 1.234567;
                static const array<array<double, 3>, 4> y
                {{
                    array<double, 3>{{0.280209,   0.909076,   0.191122}},
                    array<double, 3>{{0.650764,   0.967410,   0.044342}},
                    array<double, 3>{{0.474950,   0.838596,   0.210205}},
                    array<double, 3>{{0.421354,   0.060448,   0.267996}},
                }};
                static const array<array<double, 3>, 4> z
                {{
                    array<double, 3>{{0.345937,   1.122315,   0.235953}},
                    array<double, 3>{{0.803412,   1.194332,   0.054743}},
                    array<double, 3>{{0.586358,   1.035303,   0.259512}},
                    array<double, 3>{{0.520190,   0.074627,   0.330859}},
                }};

                array<array<double, 3>, 4> result = x * y;

                for (unsigned i(0) ; i < 4 ; ++i)
                {
                    for (unsigned j(0) ; j < 3 ; ++j)
                    {
                        TEST_CHECK_NEARLY_EQUAL(result[i][j], z[i][j], 1e-6);
                    }
                }
            }
        }
} matrix_multiplication_test;
