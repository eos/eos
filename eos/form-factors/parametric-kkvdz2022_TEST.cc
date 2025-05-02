/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2025 Danny van Dyk
 * Copyright (c) 2022      Stephan Kürten
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

Parameters
subtraction_parameters()
{
    Parameters result = Parameters::Defaults();
    result["B->gamma^*::c_1_0@KKvDZ2022"] = 0.1;
    result["B->gamma^*::c_1_1@KKvDZ2022"] = 0.2;
    result["B->gamma^*::c_1_2@KKvDZ2022"] = 0.3;
    result["B->gamma^*::c_2_0@KKvDZ2022"] = -0.1;
    result["B->gamma^*::c_2_1@KKvDZ2022"] = -0.2;
    result["B->gamma^*::c_2_2@KKvDZ2022"] = -0.3;
    result["B->gamma^*::c_3_0@KKvDZ2022"] = 0.121;
    result["B->gamma^*::c_3_1@KKvDZ2022"] = -0.2;
    result["B->gamma^*::c_3_2@KKvDZ2022"] = 0.2313;
    result["B->gamma^*::c_4_0@KKvDZ2022"] = 0.91;
    result["B->gamma^*::c_4_1@KKvDZ2022"] = -0.2;
    result["B->gamma^*::c_4_2@KKvDZ2022"] = -0.3;
    result["B->gamma^*::s_0@KKvDZ2022"] = 0.4;

    return result;
}

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
            Parameters sub_p = subtraction_parameters();

            {
                KKvDZ2022FormFactors kkvdz2022(p, Options{ });
                KKvDZ2022FormFactors kkvdz2022_sub(sub_p, Options{ {"subtracted"_ok, "on"} });

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

                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.040315453106118, 0.0               ),  kkvdz2022_sub.F_1(0.0,0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.157879067495966, 0.0824559890139291),  kkvdz2022_sub.F_1(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.157879067495966, 0.0824559890139291),  kkvdz2022_sub.F_1(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.161813624381275, 0.0835854395822285),  kkvdz2022_sub.F_1(1.2,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.143768399179354, 0.0668158087758451),  kkvdz2022_sub.F_1(1.5,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.147509153756751, 0.0677310119135683),  kkvdz2022_sub.F_1(1.5,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.683755659866785, 0.6224967477184292),  kkvdz2022_sub.F_1(0.6,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.198103842587939, 0.1328424699649706),  kkvdz2022_sub.F_1(0.9,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 0.691549948563531, 0.6311919389614886),  kkvdz2022_sub.F_1(0.6,1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-0.202648396227549, 0.1346919207079736),  kkvdz2022_sub.F_1(0.9,1.5),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.5966688305961358, 0.0                ),  kkvdz2022_sub.F_2(0.0,0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.0064639347474795, -1.0002473383097592),  kkvdz2022_sub.F_2(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.0064639347474795, -1.0002473383097592),  kkvdz2022_sub.F_2(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.0490987885176133, -1.0137950850113697),  kkvdz2022_sub.F_2(1.2,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 2.8364527078961395, -0.8105530723052206),  kkvdz2022_sub.F_2(1.5,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 2.8767700122952955, -0.8215311180858934),  kkvdz2022_sub.F_2(1.5,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.077674203759454 , -7.485189551456463 ),  kkvdz2022_sub.F_2(0.6,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.4804602984526265, -1.6108308644115616),  kkvdz2022_sub.F_2(0.9,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.17500027605132  , -7.589170522606418 ),  kkvdz2022_sub.F_2(0.6,1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.5304034457122113, -1.6330310861819934),  kkvdz2022_sub.F_2(0.9,1.5),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>( 1.817668830596136 , 0.0                ),  kkvdz2022_sub.F_3(0.0,0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.2425465832942786, -1.005030086665565 ),  kkvdz2022_sub.F_3(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.2425465832942786, -1.005030086665565 ),  kkvdz2022_sub.F_3(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.2933412747782   , -1.0211626105737832),  kkvdz2022_sub.F_3(1.2,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.0717020420752887, -0.8144282267341435),  kkvdz2022_sub.F_3(1.5,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.1197286351471765, -0.8275005464793553),  kkvdz2022_sub.F_3(1.5,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-6.929571280799387 , -7.561035286914558 ),  kkvdz2022_sub.F_3(0.6,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.737244829324993 , -1.626652973333924 ),  kkvdz2022_sub.F_3(0.9,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(-7.047636441789464 , -7.686588268074507 ),  kkvdz2022_sub.F_3(0.6,1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>( 3.7973725531421114, -1.653349634358173 ),  kkvdz2022_sub.F_3(0.9,1.5),  eps);

                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.110147091121624  , 0.0                ),  kkvdz2022_sub.F_4(0.0,0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.2794201400568068 , -0.1183084144539568),  kkvdz2022_sub.F_4(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.2794201400568068 , -0.1183084144539568),  kkvdz2022_sub.F_4(1.2,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.2859256813715785 , -0.1202591567237934),  kkvdz2022_sub.F_4(1.2,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.2591953094874802 , -0.0958683591780728),  kkvdz2022_sub.F_4(1.5,0.6),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.2653655441667946 , -0.0974490459246675),  kkvdz2022_sub.F_4(1.5,0.9),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(0.06903406359778763, -0.8970804192739344),  kkvdz2022_sub.F_4(0.6,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.339112662949608  , -0.1916659683324040),  kkvdz2022_sub.F_4(0.9,1.2),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(0.05503920070886814, -0.9123143136327135),  kkvdz2022_sub.F_4(0.6,1.5),  eps);
                TEST_CHECK_NEARLY_EQUAL( complex<double>(1.3467593608194122 , -0.19489911545919864),  kkvdz2022_sub.F_4(0.9,1.5),  eps);
            }
        }
} parametric_kkvdz2022_form_factor_test;
