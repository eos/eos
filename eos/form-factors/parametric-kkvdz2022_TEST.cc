/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Stephan Kürten
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
#include <eos/form-factors/parametric-kkvdz2022.hh>

using namespace test;
using namespace eos;

class ParametricKKvDZ2022FormFactorTest :
    public TestCase
{
    public:
        ParametricKKvDZ2022FormFactorTest() :
            TestCase("parametric_kkvdz2022_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-5;

            Parameters p = Parameters::Defaults();

            {
                KKvDZ2022FormFactors kkvdz2022(p, Options{ });

                // Tested against an IPython/Jupyter notebook implementation of Stephan Kürten
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.0713            , 0.0               ),  kkvdz2022.F_1(0.0, 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.6221950277343473, 0.6435132021736676),  kkvdz2022.F_1(0.6, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.1210071342797236, 0.0533904307432584),  kkvdz2022.F_1(0.9, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.6309148383537667, 0.6525028239395795),  kkvdz2022.F_1(0.6, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.1226949497043051, 0.0541337437915994),  kkvdz2022.F_1(0.9, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.0654480573515632, 0.0207051578485354),  kkvdz2022.F_1(1.2, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.0663460587566631, 0.0209887701778934),  kkvdz2022.F_1(1.2, 0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.0449622325887463, 0.0122004249494833),  kkvdz2022.F_1(1.5, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.0455791122008661, 0.0123675395421494),  kkvdz2022.F_1(1.5, 0.9),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.862            , 0.0               ),  kkvdz2022.F_2(0.0, 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.418157498653598, -7.735993788128181),  kkvdz2022.F_2(0.6, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.460417212010387, -0.647396133550685),  kkvdz2022.F_2(0.9, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.522039101130654, -7.843483844793416),  kkvdz2022.F_2(0.6, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.480634220081774, -0.656318569546043),  kkvdz2022.F_2(0.9, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.790537270678592, -0.251165431732555),  kkvdz2022.F_2(1.2, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.801288017768530, -0.254567344659744),  kkvdz2022.F_2(1.2, 0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.543179958726693, -0.148004286543728),  kkvdz2022.F_2(1.5, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.550565692299287, -0.150008853675771),  kkvdz2022.F_2(1.5, 0.9),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.862            ,  0.0              ),  kkvdz2022.F_3(0.0, 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.495677655406697, -7.814451669916108),  kkvdz2022.F_3(0.6, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.475015899392543, -0.653755417656518),  kkvdz2022.F_3(0.9, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.621624816124091, -7.944257202809012),  kkvdz2022.F_3(0.6, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.499383593136119, -0.664485082794080),  kkvdz2022.F_3(0.9, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.794376876741049, -0.252366435151786),  kkvdz2022.F_3(1.2, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.807203260126722, -0.256417415901600),  kkvdz2022.F_3(1.2, 0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.545816609930651, -0.148711894119271),  kkvdz2022.F_3(1.5, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.554627664985141, -0.151098878207787),  kkvdz2022.F_3(1.5, 0.9),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.1017           ,  0.0              ),  kkvdz2022.F_4(0.0, 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.895579061302733, -0.927335159257834),  kkvdz2022.F_4(0.6, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.174473803904207, -0.077031897597329),  kkvdz2022.F_4(0.9, 1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.910889380039254, -0.943085894653229),  kkvdz2022.F_4(0.6, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.177428020554160, -0.078331337250004),  kkvdz2022.F_4(0.9, 1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.093843211071694, -0.029707860390286),  kkvdz2022.F_4(1.2, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.095395895472241, -0.030197705486835),  kkvdz2022.F_4(1.2, 0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.064471059532659, -0.017505340990020),  kkvdz2022.F_4(1.5, 0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.065537628075751, -0.017793972290078),  kkvdz2022.F_4(1.5, 0.9),  eps);
            }
        }
} parametric_kkvdz2022_form_factor_test;
