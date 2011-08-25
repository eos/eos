/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <eos/utils/complex.hh>
#include <eos/utils/matrix.hh>

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

            // matrix times matrix, double
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

            // matrix times scalar, double
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

            // matrix times matrix, complex<double>
            {
                static const array<array<complex<double>, 4>, 6> x
                {{
                    array<complex<double>, 4>{{0.256065,   0.277201,   0.406745,   0.188430}},
                    array<complex<double>, 4>{{0.099332,   0.576983,   0.077084,   0.279125}},
                    array<complex<double>, 4>{{0.433378,   0.987969,   0.198432,   0.497537}},
                    array<complex<double>, 4>{{0.795675,   0.253245,   0.740202,   0.958491}},
                    array<complex<double>, 4>{{0.732521,   0.606855,   0.893422,   0.532790}},
                    array<complex<double>, 4>{{0.935207,   0.889252,   0.079741,   0.098048}},
                }};
                static const array<array<complex<double>, 3>, 4> y
                {{
                    array<complex<double>, 3>{{0.280209,   0.909076,   0.191122}},
                    array<complex<double>, 3>{{0.650764,   0.967410,   0.044342}},
                    array<complex<double>, 3>{{0.474950,   0.838596,   0.210205}},
                    array<complex<double>, 3>{{0.421354,   0.060448,   0.267996}},
                }};
                static const array<array<complex<double>, 3>, 6> z
                {{
                    array<complex<double>, 3>{{0.52472,   0.85343,   0.19723}},
                    array<complex<double>, 3>{{0.55753,   0.72999,   0.13558}},
                    array<complex<double>, 3>{{1.06826,   1.54622,   0.30169}},
                    array<complex<double>, 3>{{1.14318,   1.64699,   0.57577}},
                    array<complex<double>, 3>{{1.24900,   2.03442,   0.49750}},
                    array<complex<double>, 3>{{0.91993,   1.78324,   0.26121}},
                }};

                array<array<complex<double>, 3>, 6> result = x * y;

                for (unsigned i(0) ; i < 6 ; ++i)
                {
                    for (unsigned j(0) ; j < 3 ; ++j)
                    {
                        TEST_CHECK_NEARLY_EQUAL(result[i][j], z[i][j], 1e-5);
                    }
                }
            }

            // matrix times scalar, complex<double>
            {
                static const complex<double> x = complex<double>(1.234567, 0.3214);
                static const array<array<complex<double>, 3>, 4> y
                {{
                    array<complex<double>, 3>{{0.280209,   0.909076,   0.191122}},
                    array<complex<double>, 3>{{0.650764,   0.967410,   0.044342}},
                    array<complex<double>, 3>{{0.474950,   0.838596,   0.210205}},
                    array<complex<double>, 3>{{0.421354,   0.060448,   0.267996}},
                }};
                static const array<array<complex<double>, 3>, 4> z
                {{
                    array<complex<double>, 3>{{complex<double>(0.345936,0.090059),complex<double>(1.122315,0.292177),complex<double>(0.235953,0.061426)}},
                    array<complex<double>, 3>{{complex<double>(0.803411,0.209155),complex<double>(1.194332,0.310925),complex<double>(0.054743,0.014251)}},
                    array<complex<double>, 3>{{complex<double>(0.586358,0.152649),complex<double>(1.035303,0.269525),complex<double>(0.259512,0.067559)}},
                    array<complex<double>, 3>{{complex<double>(0.520189,0.135423),complex<double>(0.074627,0.019427),complex<double>(0.330859,0.086133)}},
                }};

                array<array<complex<double>, 3>, 4> result = x * y;

                for (unsigned i(0) ; i < 4 ; ++i)
                {
                    for (unsigned j(0) ; j < 3 ; ++j)
                    {
                        TEST_CHECK_NEARLY_EQUAL(result[i][j], z[i][j], 1e-6);
                    }
                }
            }

            // vector * matrix
            {
                typedef array<array<double, 3>, 3> Matrix;
                typedef array<double, 3> Vector;

                const Matrix A
                {{
                    {{1, 2, 3}},
                    {{4, 5, 6}},
                    {{7, 8, 9}}
                }};

                const Vector x {{1/2., 1/3., 1/4.}};

                // result
                const Vector true_result {{43/12., 14/3., 23/4.}};
                Vector y = x ^ A;

                for (unsigned i = 0 ; i < 3 ; ++i)
                    TEST_CHECK_RELATIVE_ERROR(y[i], true_result[i], 1e-15);
            }

            // scalar product
            {
                typedef array<double, 3> Vector;

                const Vector x {{1/2., 1/3., 1/4.}};
                const Vector y {{43/12., 14/3., 23/4.}};

                TEST_CHECK_RELATIVE_ERROR(dot(x,y), 4.7847222222222222, 1e-15);
            }

            // vector - vector
            {
                typedef array<double, 3> Vector;

                const Vector x {{1 / 2., 1 / 3., 1 / 4.}};
                const Vector y {{43 / 12., 14 / 3., 23 / 4.}};

                Vector true_result {{-37 / 12., -52 / 12., -66 / 12.}};
                Vector result = x - y;

                for (unsigned i = 0 ; i < 3 ; ++i)
                {
                    TEST_CHECK_RELATIVE_ERROR(result[i], true_result[i], 1e-15);
                }
            }
        }
} matrix_multiplication_test;
