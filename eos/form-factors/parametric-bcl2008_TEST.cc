/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2014, 2015, 2018 Danny van Dyk
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
#include <eos/form-factors/parametric-bcl2008-impl.hh>

using namespace test;
using namespace eos;

class BCL2008FormFactorsTest :
    public TestCase
{
    public:
        BCL2008FormFactorsTest() :
            TestCase("bcl2008_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::BCL2008", p, Options{ });

                p["B->pi::f_+(0)@BCL2008"] = 1.0;
                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.39360, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19462, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.48329, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.95593, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.87063, eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.16483, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.40109, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.77364, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 2.47348, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.17917, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.45663, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.94890, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.06978, eps);

                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.14259, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.33718, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.62004, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.07054, eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.16740, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.41489, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.82327, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.64024, eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.18391, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.45892, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.91416, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.80533, eps);
            }
        }
} bcl2008_form_factors_test;

class BCL2008K4FormFactorsTest :
    public TestCase
{
    public:
        BCL2008K4FormFactorsTest() :
            TestCase("bcl2008_K4_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::BCL2008-4", p, Options{ });

                p["B->pi::f_+(0)@BCL2008"] = 1.0;
                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.39360, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19462, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.48329, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.95593, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.87063, eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.16026, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.39059, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.75519, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 2.44228, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.14259, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.33718, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.62004, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.07054, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.18833, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.47763, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.98580, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.13217, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.16740, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.41489, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.82327, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.64024, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.20035, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.51330, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.06777, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.30001, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.18391, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.45892, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.91416, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.80533, eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19087, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.47546, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.94362, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.85213, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.39360, eps);

                p["B->pi::b_T^1@BCL2008"]  = 1.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.16026, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.39059, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 1.75519, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 2.44228, eps);

                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 1.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.18833, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.47763, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 1.98580, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.13217, eps);

                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.20035, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.51330, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.06777, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.30001, eps);
            }
        }
} bcl2008_k4_form_factors_k4_test;

class BCL2008K5FormFactorsTest :
    public TestCase
{
    public:
        BCL2008K5FormFactorsTest() :
            TestCase("bcl2008_K5_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::BCL2008-5", p, Options{ });

                p["B->pi::f_+(0)@BCL2008"] = 1.0;
                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.39360, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19462, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.48329, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.95593, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.87063, eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.16146, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.39313, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.75929, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 2.44899, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.14259, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.33718, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.62004, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.07054, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.18592, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.47256, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.97759, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.11876, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.16740, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.41489, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.82327, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.64024, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 1.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.20396, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.52090, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.08008, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.32013, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.18391, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.45892, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.91416, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.80533, eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 1.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 1.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.20927, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.53465, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.10670, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.36677, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19087, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.47546, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.94362, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.85213, eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.19338, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.48090, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.95239, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.86539, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.39360, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 1.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.16146, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.39313, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 1.75929, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 2.44899, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 1.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.18592, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.47256, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 1.97759, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.11876, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 1.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.20396, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.52090, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.08008, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.32013, eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.20927, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.53465, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.10670, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 3.36677, eps);
            }
        }
} bcl2008_k5_form_factors_k5_test;
