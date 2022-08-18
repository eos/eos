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

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19462, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.48329, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.95593, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.87063, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16483, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.40109, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.77364, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.47348, ff->f_p(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.17917, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.45663, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.94890, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.06978, ff->f_p(20.0), eps);

                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.14259, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.33718, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.62004, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.07054, ff->f_0(20.0), eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16740, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.41489, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.82327, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.64024, ff->f_0(20.0), eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18391, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.45892, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.91416, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.80533, ff->f_0(20.0), eps);
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

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19462, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.48329, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.95593, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.87063, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16026, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.39059, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.75519, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.44228, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.14259, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.33718, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.62004, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.07054, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18833, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47763, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.98580, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.13217, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16740, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.41489, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.82327, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.64024, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20035, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.51330, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.06777, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.30001, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18391, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.45892, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.91416, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.80533, ff->f_0(20.0), eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19087, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47546, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.94362, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.85213, ff->f_0(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_t(20.0), eps);

                p["B->pi::b_T^1@BCL2008"]  = 1.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16026, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.39059, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.75519, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.44228, ff->f_t(20.0), eps);

                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 1.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18833, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47763, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.98580, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.13217, ff->f_t(20.0), eps);

                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20035, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.51330, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.06777, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.30001, ff->f_t(20.0), eps);
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

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19462, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.48329, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.95593, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.87063, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 1.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 1.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16146, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.39313, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.75929, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.44899, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.14259, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.33718, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.62004, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.07054, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 1.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 1.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18592, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47256, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.97759, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.11876, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16740, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.41489, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.82327, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.64024, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 1.0;
                p["B->pi::b_+^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 1.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20396, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.52090, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.08008, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.32013, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18391, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.45892, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.91416, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.80533, ff->f_0(20.0), eps);

                p["B->pi::b_+^1@BCL2008"]  = 0.0;
                p["B->pi::b_+^2@BCL2008"]  = 0.0;
                p["B->pi::b_+^3@BCL2008"]  = 0.0;
                p["B->pi::b_+^4@BCL2008"]  = 1.0;
                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 1.0;
                p["B->pi::b_0^5@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20927, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.53465, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.10670, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.36677, ff->f_p(20.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19087, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47546, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.94362, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.85213, ff->f_0(20.0), eps);

                p["B->pi::b_0^1@BCL2008"]  = 0.0;
                p["B->pi::b_0^2@BCL2008"]  = 0.0;
                p["B->pi::b_0^3@BCL2008"]  = 0.0;
                p["B->pi::b_0^4@BCL2008"]  = 0.0;
                p["B->pi::b_0^5@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.19338, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.48090, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.95239, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.86539, ff->f_0(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_t(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 1.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16146, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.39313, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.75929, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.44899, ff->f_t(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 1.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.18592, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.47256, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.97759, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.11876, ff->f_t(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 1.0;
                p["B->pi::b_T^4@BCL2008"]  = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20396, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.52090, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.08008, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.32013, ff->f_t(20.0), eps);

                p["B->pi::f_T(0)@BCL2008"] = 1.0;
                p["B->pi::b_T^1@BCL2008"]  = 0.0;
                p["B->pi::b_T^2@BCL2008"]  = 0.0;
                p["B->pi::b_T^3@BCL2008"]  = 0.0;
                p["B->pi::b_T^4@BCL2008"]  = 1.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.20927, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.53465, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.10670, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.36677, ff->f_t(20.0), eps);
            }
        }
} bcl2008_k5_form_factors_k5_test;
