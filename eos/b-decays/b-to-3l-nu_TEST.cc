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
#include <eos/b-decays/b-to-3l-nu.hh>

using namespace test;
using namespace eos;

class BToThreeLeptonsNeutrinoTest :
    public TestCase
{
    public:
        BToThreeLeptonsNeutrinoTest() :
            TestCase("b_to_3l_nu_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();

            {

                BToThreeLeptonsNeutrino three_lnu(p, Options{ { "l"_ok, "mu" }, { "lprime"_ok, "e" } });

                // Tested against an IPython/Jupyter notebook implementation of Stephan Kürten
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.0,  0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.0,  1.5), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 1.5,  0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 1.2,  0.6), 4.39594668126115e-11  , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 1.2,  0.9), 4.470741659137669e-11 , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 1.5,  0.6), 1.950546024128868e-11 , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 1.5,  0.9), 1.9810430065283566e-11, 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.6,  1.2), 7.984365465773253e-09 , 1e-22);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.9,  1.2), 1.6702367962034642e-10, 1e-23);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.6,  1.5), 8.177374478468938e-09 , 1e-22);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.9,  1.5), 1.7009605852663074e-10, 1e-23);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio( 0.0, 50.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_branching_ratio(50.0,  1.5), 0.0                   , 1e-24);
                /*
                 * differential decay width of 5 kinematic variables
                 * quintuple_differential_decay_width(q2, k2, z_gamma, z_w, phi)
                 * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
                 * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
                 * z_gamma is the angle between the negatively charged lepton l' and the negative z-axis
                 * z_w is the angle between the charged lepton l and the positive z-axis
                 * phi is the angle between the q2 plane and the k2 plane
                 */
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.0,  0.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.0,  1.5, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 1.5,  0.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 0.0, 0.0, 0.0), 1.4012775797795022e-11, 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 1.0,-1.0, 3.1), 5.1287437263386674e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 1.0,-1.0,-3.1), 5.1287437263386674e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5,-1.0,-1.0,-3.1), 5.1287437263386674e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 1.0, 1.0, 3.1), 3.085696796294097e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 0.6, 0.4, 3.1), 1.1160400109072714e-11, 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 0.4, 0.6, 0.3), 6.643678742005261e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.9,  1.5, 0.4, 0.6,-0.3), 6.644934752365039e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.6,  1.2, 0.4, 0.6,-0.3), 3.1035668687398013e-10, 1e-23);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 1.2,  0.6, 0.4, 0.6,-0.3), 1.9082580462960103e-12, 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio( 0.0, 50.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.quintuple_differential_branching_ratio(50.0,  1.5, 0.0, 0.0, 0.0), 0.0                   , 1e-24);

                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 0.0,   0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 0.0,   1.5), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 1.5,   0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 1.0,   4.7), -0.19978703404850154, 1e-13);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 0.6,   3.2), -0.1900337624486565 , 1e-13);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 0.9,   3.2), -0.15573648794445794, 1e-13);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 1.5,   4.7), -0.1703444419281766 , 1e-13);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 1.2,   4.7), -0.18574472584759466, 1e-13);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry( 0.0,  50.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.double_differential_forward_backward_asymmetry(50.0,   1.5), 0.0                   , 1e-24);

                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_branching_ratio(1.0, 4.0, 1.0, 4.0), 1.11942599097845e-10  , 1e-14);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_branching_ratio(4.0, 5.0, 3.0, 5.0), 1.8054715753057808e-12, 1e-16);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_branching_ratio(1.2, 3.0, 0.2, 5.0), 1.0192934595042032e-10, 1e-14);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_branching_ratio(1.2, 3.5, 2.2, 5.0), 6.514881129640035e-11 , 1e-15);

                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_forward_backward_asymmetry(1.0, 4.0, 1.0, 4.0), -0.10720983486749855, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_forward_backward_asymmetry(4.0, 5.0, 3.0, 5.0), -0.09841640103003944, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_forward_backward_asymmetry(1.2, 3.0, 0.2, 5.0), -0.10745200783863061, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu.integrated_forward_backward_asymmetry(1.2, 3.5, 2.2, 5.0), -0.1334186883151141 , 1e-5);

                BToThreeLeptonsNeutrino three_lnu_with_tau(p, Options{ { "l"_ok, "tau" }, { "lprime"_ok, "mu" } });

                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.0,  0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.0,  1.5), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 1.5,  0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 1.2,  3.2), 2.2656835615631178e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 1.5,  3.2), 9.77022421527314e-15  , 1e-28);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 1.2,  3.5), 1.1407883572580513e-12, 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 1.5,  3.5), 4.900016933078733e-13 , 1e-26);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.6,  4.1), 1.1431202742303967e-09, 1e-22);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.9,  4.1), 2.1967778054652107e-11, 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.6, 15.0), 9.798358180136963e-09 , 1e-22);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.9, 15.0), 1.4812669954450616e-10, 1e-23);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio( 0.0, 50.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_branching_ratio(50.0,  1.5), 0.0                   , 1e-24);

                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.0,  0.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.0,  1.5, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 1.5,  0.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 0.0, 0.0, 0.0), 2.031669886855895e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 1.0,-1.0, 3.1), 1.7424121078854753e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 1.0,-1.0,-3.1), 1.7424121078854753e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5,-1.0,-1.0,-3.1), 1.7424121078854753e-14, 1e-27);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 1.0, 1.0, 3.1), 7.81354285945554e-13  , 1e-26);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 0.6, 0.4, 3.1), 1.952956580972857e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 0.4, 0.6, 0.3), 2.5081122681564525e-12, 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.9,  4.5, 0.4, 0.6,-0.3), 2.5075950200134066e-12, 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.6,  4.2, 0.4, 0.6,-0.3), 9.292842351120294e-11 , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 1.2, 14.6, 0.4, 0.6,-0.3), 1.500147195568813e-12 , 1e-25);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio( 0.0, 50.0, 0.0, 0.0, 0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.quintuple_differential_branching_ratio(50.0,  1.5, 0.0, 0.0, 0.0), 0.0                   , 1e-24);

                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 0.0,   0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 0.0,   1.5), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 1.5,   0.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 1.0, 4.7), -0.43890105637843724, 1e-10);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 0.6, 3.2), -0.49813547368786776, 1e-10);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 0.9, 3.2), -0.4958929059233193 , 1e-10);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 1.5, 4.7), -0.41899237708674814, 1e-10);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 1.2, 4.7), -0.4307210623669162 , 1e-10);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry( 0.0,50.0), 0.0                   , 1e-24);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.double_differential_forward_backward_asymmetry(50.0, 1.5), 0.0                   , 1e-24);

                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_branching_ratio(1.0, 4.0, 3.5, 4.0), 1.0814431951945232e-12, 1e-16);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_branching_ratio(4.0, 5.0, 3.2, 5.0), 1.2892379916914902e-13, 1e-17);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_branching_ratio(1.2, 3.0, 3.5, 5.0), 4.327260506854218e-12 , 1e-16);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_branching_ratio(1.2, 3.5, 3.2, 5.0), 4.568873862504263e-12 , 1e-16);

                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_forward_backward_asymmetry(1.0, 4.0, 3.5, 4.0), -0.4580920453975816 , 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_forward_backward_asymmetry(4.0, 5.0, 3.2, 5.0), -0.3173962670686166 , 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_forward_backward_asymmetry(1.2, 3.0, 3.5, 5.0), -0.42333712608602736, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(three_lnu_with_tau.integrated_forward_backward_asymmetry(1.2, 3.5, 3.2, 5.0), -0.4218080496815521 , 1e-5);
            }
        }
} b_to_3l_nu_test;
