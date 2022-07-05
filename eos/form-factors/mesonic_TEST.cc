/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2014, 2015, 2018 Danny van Dyk
 * Copyright (c) 2015 Frederik Beaujean
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
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
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic-impl.hh>

#include <vector>

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

class BSZ2015FormFactorsTest :
    public TestCase
{
    public:
        BSZ2015FormFactorsTest() :
            TestCase("bsz2015_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                // test case created by using the known relations among the BCL2008 parameters
                // for the highest power.
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::BSZ2015", p, Options{ });

                p["B->pi::alpha^f+_0@BSZ2015"] = 1.0;
                p["B->pi::alpha^f+_1@BSZ2015"] = 0.0;
                p["B->pi::alpha^f+_2@BSZ2015"] = 0.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);

                p["B->pi::alpha^f+_1@BSZ2015"] = 1.0;
                p["B->pi::alpha^f+_2@BSZ2015"] = 2.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16581, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.42261, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.88375, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.97499, ff->f_p(20.0), eps);

                p["B->pi::alpha^f0_1@BSZ2015"] = 1.5;
                p["B->pi::alpha^f0_2@BSZ2015"] = 1.5;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.11998, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.28572, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.53862, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.00499, ff->f_0(20.0), eps);

                p["B->pi::alpha^fT_0@BSZ2015"] =  1.0;
                p["B->pi::alpha^fT_1@BSZ2015"] = -1.0;
                p["B->pi::alpha^fT_2@BSZ2015"] =  2.5;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.27271, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.73442, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.64425, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(4.99850, ff->f_t(20.0), eps);
            }
        }
} bsz2015_form_factors_test;

class BToKstarKMPW2010FormFactorsTest:
    public TestCase
{
    public:
        BToKstarKMPW2010FormFactorsTest() :
            TestCase("b_to_kstar_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5e-6;

            // central values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.36; p["B->K^*::b^V_1@KMPW2010"]  =  -4.8;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.29; p["B->K^*::b^A0_1@KMPW2010"] = -18.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.25; p["B->K^*::b^A1_1@KMPW2010"] =  +0.34;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.23; p["B->K^*::b^A2_1@KMPW2010"] =  -0.85;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.31; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.31; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.22; p["B->K^*::b^T3_1@KMPW2010"] = -10.3;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.421964, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.497285, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.590298, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.707366, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.858208, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.058404, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.334545, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.735915, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.410875, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.560285, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.747894, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.988041, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.302880, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.728452, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.327350, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(3.218300, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.266631, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.285779, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.308075, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334383, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.365917, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.404440, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.452618, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514681, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.250198, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.273759, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.301563, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334822, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.375252, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.425358, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.488950, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572108, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.362235, ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.425698, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.504029, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.602573, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.729487, ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.897853, ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.12999 , ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.46727 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.35026 , ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.397838, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.45472 , ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.523652, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.608536, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.715103, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.852101, ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.03355 , ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.276523, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.344214, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.426193, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.526794, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.6522  , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.811519, ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.01871 , ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.29621 , ff->t_3(18.4), eps);
            }

            // raised values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.59; p["B->K^*::b^V_1@KMPW2010"]  =  -4.0;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.39; p["B->K^*::b^A0_1@KMPW2010"] = -16.9;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.41; p["B->K^*::b^A1_1@KMPW2010"] =  +1.2;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.42; p["B->K^*::b^A2_1@KMPW2010"] =  +2.03;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.49; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 + 0.81;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.49; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 + 2.1;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.39; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 + 2.5;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.683006, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.795816, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.934831, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.109443, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.333989, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.631445, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.041006, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.635294, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.543363, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.732820, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.970583, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.274770, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.673378, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(2.211932, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.969510, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(4.096084, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.430965, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.454709, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.481882, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.513361, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.550366, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.594652, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.648841, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.717050, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.435238, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.451990, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.470539, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.491252, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514616, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.541305, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572274, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.608954, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.56538 , ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.656751, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.769278, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.910534, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.09208 , ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.33244 , ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.66321 , ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.14292 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.535223, ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.588078, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.650576, ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.725486, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.816729, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.930041, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.07414 , ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.26297 , ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.472753, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.571573, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.690926, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.837005, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.01864 , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.24883 , ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.54748 , ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.94657 , ff->t_3(18.4), eps);            }

            // lowered values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.24; p["B->K^*::b^V_1@KMPW2010"]  =  -5.2;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.22; p["B->K^*::b^A0_1@KMPW2010"] = -21.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.15; p["B->K^*::b^A1_1@KMPW2010"] =  -0.46;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.13; p["B->K^*::b^A2_1@KMPW2010"] =  -2.2;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.21; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 - 0.41;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.21; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 - 2.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.12; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 - 3.1;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.2830474, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3354240, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4001625, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4817170, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.5868880, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.7265850, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.9194250, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.1999232, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.3236650, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4519470, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.6131980, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.8198090, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.0909360, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.4577380, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.9743500, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.7434500, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.1621258, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1762210, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1927954, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2125507, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2364774, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2660210, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3033756, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3520390, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.1445568, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1616854, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1820760, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2066815, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2368550, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2745790, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3228774, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3865910, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.246944, ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.291875, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.347386, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.417288, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.507397, ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.627041, ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.792141, ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.03221 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.24554 , ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.287805, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.338645, ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.400627, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.477404, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.574348, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.699678, ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.866585, ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.157488, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.202489, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.257115, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.324297, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.40822 , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.515053, ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.654254, ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.841037, ff->t_3(18.4), eps);
            }
        }
} b_to_kstar_kmpw2010_form_factors_test;

class BToKKMPW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKKMPW2010FormFactorsTest() :
            TestCase("b_to_k_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K::KMPW2010", p, Options{ });
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.3844144474375, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4381960494587, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5043916755865, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5874824689246, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6943252577904, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8359400978914, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.0311821367700, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.3150955485260, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(1.7607545970800, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.4417893713152, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5045402114062, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5818230155226, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6788902047561, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8037805019019, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9694161202369, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.1979110863960, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.5303740143290, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.0525389704830, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3691335300041, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4001522415666, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4333067770539, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4689025169045, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5073172894003, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5490271105427, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5946449022903, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6449812193912, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7011445924499, ff->f_0(20.7), eps);
        }
} b_to_k_kmpw2010_form_factors_test;

class BToBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToBFW2010FormFactorsTest() :
            TestCase("b_to_k_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K::BFW2010", p, Options{ });
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.3285329106071, ff->f_p( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3765270351772, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4352582092738, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5083556849771, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6012033334785, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7221134890713, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8845704673796, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.1118693561113, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.4477383401655, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(1.9842848427703, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3529456128368, ff->f_t( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4057232066867, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4703717236260, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5509131505349, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6533164641168, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7868011822009, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9663293753277, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.2177594231842, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.5896500926049, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.1843219169364, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3208496768629, ff->f_0( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3456358228669, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3724574045224, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4016438048517, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4336098586343, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4688881948925, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5081785681054, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5524262352940, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6029528628018, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6616893406433, ff->f_0(20.7), eps);
        }
} b_to_k_bfw2010_form_factors_test;

class BToDstarBSZ2015FormFactorsTest :
public TestCase
{
public:
    BToDstarBSZ2015FormFactorsTest() :
    TestCase("b_to_Dstar_bsz2015_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 5.1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->D^*::BSZ2015", p, Options{ });
        TEST_CHECK(ff.get() != nullptr);

        /*
         * for the B->D^* case no SSE parameters are available yet, so for the moment we set them to zero and once the GKvD2018 calculation is fitted we can then update the test case later on.
         */
        p["B->D^*::alpha^A0_0@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A0_1@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A0_2@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A1_0@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A1_1@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A1_2@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^A12_1@BSZ2015"] = 0.0;
        p["B->D^*::alpha^A12_2@BSZ2015"] = 0.0;
        p["B->D^*::alpha^V_0@BSZ2015"  ] = 0.0;
        p["B->D^*::alpha^V_1@BSZ2015"  ] = 0.0;
        p["B->D^*::alpha^V_2@BSZ2015"  ] = 0.0;
        p["B->D^*::alpha^T1_0@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^T1_1@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^T1_2@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^T2_1@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^T2_2@BSZ2015" ] = 0.0;
        p["B->D^*::alpha^T23_0@BSZ2015"] = 0.0;
        p["B->D^*::alpha^T23_1@BSZ2015"] = 0.0;
        p["B->D^*::alpha^T23_2@BSZ2015"] = 0.0;

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(6.1), eps);
    }
} b_to_Dstar_bsz2015_form_factors_test;

class BToKstarBSZ2015FormFactorsTest :
    public TestCase
{
    public:
        BToKstarBSZ2015FormFactorsTest() :
            TestCase("b_to_kstar_bsz2015_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5.1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::BSZ2015", p, Options{ });
            TEST_CHECK(ff.get() != nullptr);

            /* compare with values from David Straub,
             * use his values of the parameters
             */
            p["B->K^*::alpha^A0_0@BSZ2015" ] = +0.39;
            p["B->K^*::alpha^A0_1@BSZ2015" ] = -1.15;
            p["B->K^*::alpha^A0_2@BSZ2015" ] = +2.08;
            p["B->K^*::alpha^A1_0@BSZ2015" ] = +0.29;
            p["B->K^*::alpha^A1_1@BSZ2015" ] = +0.31;
            p["B->K^*::alpha^A1_2@BSZ2015" ] = +0.72;
            p["B->K^*::alpha^A12_1@BSZ2015"] = +0.57;
            p["B->K^*::alpha^A12_2@BSZ2015"] = +0.14;
            p["B->K^*::alpha^V_0@BSZ2015"  ] = +0.37;
            p["B->K^*::alpha^V_1@BSZ2015"  ] = -1.08;
            p["B->K^*::alpha^V_2@BSZ2015"  ] = +2.47;
            p["B->K^*::alpha^T1_0@BSZ2015" ] = +0.31;
            p["B->K^*::alpha^T1_1@BSZ2015" ] = -0.96;
            p["B->K^*::alpha^T1_2@BSZ2015" ] = +2.01;
            p["B->K^*::alpha^T2_1@BSZ2015" ] = +0.42;
            p["B->K^*::alpha^T2_2@BSZ2015" ] = +2.02;
            p["B->K^*::alpha^T23_0@BSZ2015"] = +0.79;
            p["B->K^*::alpha^T23_1@BSZ2015"] = +1.26;
            p["B->K^*::alpha^T23_2@BSZ2015"] = +1.96;

            TEST_CHECK_NEARLY_EQUAL(0.393136, ff->a_0(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.440394, ff->a_0(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.496878, ff->a_0(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.565342, ff->a_0(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.289606, ff->a_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3039,   ff->a_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.319847, ff->a_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.337861, ff->a_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.248569, ff->a_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.272718, ff->a_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.300677, ff->a_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.333431, ff->a_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.367312, ff->v(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.411249, ff->v(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.463812, ff->v(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.527595, ff->v(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3094,   ff->t_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.346962, ff->t_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.391946, ff->t_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.446575, ff->t_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.308387, ff->t_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.322844, ff->t_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.339239, ff->t_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.358183, ff->t_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.184952, ff->t_3(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.200925, ff->t_3(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.219004, ff->t_3(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.239587, ff->t_3(6.1), eps);
        }
} b_to_kstar_bsz2015_form_factors_test;

class BToRhoBSZ2015FormFactorsTest :
public TestCase
{
public:
    BToRhoBSZ2015FormFactorsTest() :
    TestCase("b_to_rho_bsz2015_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 5.1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->rho::BSZ2015", p, Options{ });
        TEST_CHECK(ff.get() != nullptr);

        /*
         * use David Straub's values of the SSE parameters from LCSR only
         */
        p["B->rho::alpha^A0_0@BSZ2015" ] = +0.36;
        p["B->rho::alpha^A0_1@BSZ2015" ] = -0.83;
        p["B->rho::alpha^A0_2@BSZ2015" ] = +1.33;
        p["B->rho::alpha^A1_0@BSZ2015" ] = +0.26;
        p["B->rho::alpha^A1_1@BSZ2015" ] = +0.39;
        p["B->rho::alpha^A1_2@BSZ2015" ] = +0.16;
        p["B->rho::alpha^A12_1@BSZ2015"] = +0.76;
        p["B->rho::alpha^A12_2@BSZ2015"] = +0.46;
        p["B->rho::alpha^V_0@BSZ2015"  ] = +0.33;
        p["B->rho::alpha^V_1@BSZ2015"  ] = -0.86;
        p["B->rho::alpha^V_2@BSZ2015"  ] = +1.80;
        p["B->rho::alpha^T1_0@BSZ2015" ] = +0.27;
        p["B->rho::alpha^T1_1@BSZ2015" ] = -0.74;
        p["B->rho::alpha^T1_2@BSZ2015" ] = +1.45;
        p["B->rho::alpha^T2_1@BSZ2015" ] = +0.47;
        p["B->rho::alpha^T2_2@BSZ2015" ] = +0.58;
        p["B->rho::alpha^T23_0@BSZ2015"] = +0.75;
        p["B->rho::alpha^T23_1@BSZ2015"] = +1.90;
        p["B->rho::alpha^T23_2@BSZ2015"] = +2.93;

        TEST_CHECK_NEARLY_EQUAL(0.36186,  ff->a_0(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.402772, ff->a_0(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.4521,   ff->a_0(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.512422, ff->a_0(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.260532, ff->a_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.271749, ff->a_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.284225, ff->a_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.29821,  ff->a_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.226525, ff->a_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.247972, ff->a_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.272771, ff->a_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.301645, ff->a_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.331752, ff->v(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.370391, ff->v(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.417199, ff->v(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.474694, ff->v(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.271458, ff->t_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.303616, ff->t_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.342573, ff->t_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.390422, ff->t_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.270508, ff->t_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.28128,  ff->t_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.29338,  ff->t_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.307107, ff->t_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.179298, ff->t_3(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.196067, ff->t_3(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.215543, ff->t_3(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.238325, ff->t_3(6.1), eps);
    }
} b_to_rho_bsz2015_form_factors_test;


class BToPiPiFvDV2018FormFactorsTest :
    public TestCase
{
    public:
        BToPiPiFvDV2018FormFactorsTest() :
        TestCase("b_to_pi_pi_fvdv2018_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5.1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToPP>> ff = FormFactorFactory<PToPP>::create("B->pipi::FvDV2018", p, Options{ });
            TEST_CHECK(ff.get() != nullptr);

            // time
            p["B->pipi::a^Ftime_0_0@FvDV2018"] = -0.36816929771;
            p["B->pipi::a^Ftime_0_1@FvDV2018"] =  2.24240299735;
            p["B->pipi::a^Ftime_0_2@FvDV2018"] = -4.24917695998;
            p["B->pipi::a^Ftime_0_3@FvDV2018"] =  2.35541703192;
            p["B->pipi::a^Ftime_1_0@FvDV2018"] =  5.24115374042;
            p["B->pipi::a^Ftime_1_1@FvDV2018"] = -22.0437225418;
            p["B->pipi::a^Ftime_1_2@FvDV2018"] =  20.2542056471;
            p["B->pipi::b^Ftime_0_0@FvDV2018"] = -0.78770063503;
            p["B->pipi::b^Ftime_0_1@FvDV2018"] =  10.2021697105;
            p["B->pipi::b^Ftime_0_2@FvDV2018"] = -36.9473507702;
            p["B->pipi::b^Ftime_0_3@FvDV2018"] =  41.0909785695;
            p["B->pipi::b^Ftime_1_0@FvDV2018"] =  12.1068819568;
            p["B->pipi::b^Ftime_1_1@FvDV2018"] = -86.2325621190;
            p["B->pipi::b^Ftime_1_2@FvDV2018"] =  147.922534873;
            p["B->pipi::c^Ftime_0_0@FvDV2018"] =  1.56132165878;
            p["B->pipi::c^Ftime_0_1@FvDV2018"] = -15.8641613048;
            p["B->pipi::c^Ftime_0_2@FvDV2018"] =  51.0904259777;
            p["B->pipi::c^Ftime_0_3@FvDV2018"] = -52.9767513762;
            p["B->pipi::c^Ftime_1_0@FvDV2018"] = -23.6236882937;
            p["B->pipi::c^Ftime_1_1@FvDV2018"] =  141.591407772;
            p["B->pipi::c^Ftime_1_2@FvDV2018"] = -209.769464670;
            // long
            p["B->pipi::a^Flong_0_0@FvDV2018"] = -0.16139192740;
            p["B->pipi::a^Flong_0_1@FvDV2018"] =  1.09519092754;
            p["B->pipi::a^Flong_0_2@FvDV2018"] = -2.41027721486;
            p["B->pipi::a^Flong_0_3@FvDV2018"] =  1.70033143492;
            p["B->pipi::a^Flong_1_0@FvDV2018"] =  1.99824594658;
            p["B->pipi::a^Flong_1_1@FvDV2018"] = -8.62878627878;
            p["B->pipi::a^Flong_1_2@FvDV2018"] =  9.35509122572;
            p["B->pipi::b^Flong_0_0@FvDV2018"] = -0.44278949412;
            p["B->pipi::b^Flong_0_1@FvDV2018"] =  4.74822500206;
            p["B->pipi::b^Flong_0_2@FvDV2018"] = -15.9081401854;
            p["B->pipi::b^Flong_0_3@FvDV2018"] =  17.0309192822;
            p["B->pipi::b^Flong_1_0@FvDV2018"] = -5.92715117703;
            p["B->pipi::b^Flong_1_1@FvDV2018"] =  32.8787517243;
            p["B->pipi::b^Flong_1_2@FvDV2018"] = -44.5974162489;
            p["B->pipi::c^Flong_0_0@FvDV2018"] =  0.77102061709;
            p["B->pipi::c^Flong_0_1@FvDV2018"] = -6.98963227925;
            p["B->pipi::c^Flong_0_2@FvDV2018"] =  21.1956324818;
            p["B->pipi::c^Flong_0_3@FvDV2018"] = -21.2213423789;
            p["B->pipi::c^Flong_1_0@FvDV2018"] =  4.11524737745;
            p["B->pipi::c^Flong_1_1@FvDV2018"] = -23.9234305792;
            p["B->pipi::c^Flong_1_2@FvDV2018"] =  32.8488722373;
            // para
            p["B->pipi::a^Fpara_0_0@FvDV2018"] = -0.73839711162;
            p["B->pipi::a^Fpara_0_1@FvDV2018"] =  5.33340671033;
            p["B->pipi::a^Fpara_0_2@FvDV2018"] = -13.0548108152;
            p["B->pipi::a^Fpara_0_3@FvDV2018"] =  10.7823778312;
            p["B->pipi::a^Fpara_1_0@FvDV2018"] =  4.90990913468;
            p["B->pipi::a^Fpara_1_1@FvDV2018"] =  3.07689657262;
            p["B->pipi::a^Fpara_1_2@FvDV2018"] = -51.8342261186;
            p["B->pipi::b^Fpara_0_0@FvDV2018"] = -0.77364251324;
            p["B->pipi::b^Fpara_0_1@FvDV2018"] =  8.65261164441;
            p["B->pipi::b^Fpara_0_2@FvDV2018"] = -29.0610911551;
            p["B->pipi::b^Fpara_0_3@FvDV2018"] =  30.4513866623;
            p["B->pipi::b^Fpara_1_0@FvDV2018"] =  16.8847387994;
            p["B->pipi::b^Fpara_1_1@FvDV2018"] = -203.754603759;
            p["B->pipi::b^Fpara_1_2@FvDV2018"] =  456.700391055;
            p["B->pipi::c^Fpara_0_0@FvDV2018"] =  1.91667289560;
            p["B->pipi::c^Fpara_0_1@FvDV2018"] = -16.8682304272;
            p["B->pipi::c^Fpara_0_2@FvDV2018"] =  48.4938486882;
            p["B->pipi::c^Fpara_0_3@FvDV2018"] = -45.0828744220;
            p["B->pipi::c^Fpara_1_0@FvDV2018"] = -32.1123316945;
            p["B->pipi::c^Fpara_1_1@FvDV2018"] =  273.753395427;
            p["B->pipi::c^Fpara_1_2@FvDV2018"] = -531.002613990;
            // perp
            p["B->pipi::a^Fperp_0_0@FvDV2018"] = -2.05535479854;
            p["B->pipi::a^Fperp_0_1@FvDV2018"] =  16.9988604895;
            p["B->pipi::a^Fperp_0_2@FvDV2018"] = -46.6407789359;
            p["B->pipi::a^Fperp_0_3@FvDV2018"] =  42.5157927594;
            p["B->pipi::a^Fperp_1_0@FvDV2018"] =  25.4079749879;
            p["B->pipi::a^Fperp_1_1@FvDV2018"] = -139.728111016;
            p["B->pipi::a^Fperp_1_2@FvDV2018"] =  193.460283599;
            p["B->pipi::b^Fperp_0_0@FvDV2018"] =  8.23644248367;
            p["B->pipi::b^Fperp_0_1@FvDV2018"] = -67.3709309476;
            p["B->pipi::b^Fperp_0_2@FvDV2018"] =  183.591398737;
            p["B->pipi::b^Fperp_0_3@FvDV2018"] = -166.684702751;
            p["B->pipi::b^Fperp_1_0@FvDV2018"] = -112.935412057;
            p["B->pipi::b^Fperp_1_1@FvDV2018"] =  624.097066875;
            p["B->pipi::b^Fperp_1_2@FvDV2018"] = -867.134281277;
            p["B->pipi::c^Fperp_0_0@FvDV2018"] = -6.96163577513;
            p["B->pipi::c^Fperp_0_1@FvDV2018"] =  56.0519566307;
            p["B->pipi::c^Fperp_0_2@FvDV2018"] = -151.061706725;
            p["B->pipi::c^Fperp_0_3@FvDV2018"] =  136.078387363;
            p["B->pipi::c^Fperp_1_0@FvDV2018"] =  106.107678887;
            p["B->pipi::c^Fperp_1_1@FvDV2018"] = -588.669971295;
            p["B->pipi::c^Fperp_1_2@FvDV2018"] =  819.519985258;

            // time
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_time(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.068642, imag(ff->f_time(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_time(0.05, 16.0, -0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.433295, imag(ff->f_time(0.05, 16.0, -0.5)), eps);
            // long
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_long(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.032153, imag(ff->f_long(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_long(0.05, 16.0, -0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.413832, imag(ff->f_long(0.05, 16.0, -0.5)), eps);
            // para
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_para(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.021664, imag(ff->f_para(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_para(0.05, 16.0, -0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.013560, imag(ff->f_para(0.05, 16.0, -0.5)), eps);
            // perp
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_perp(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.001048, imag(ff->f_perp(0.60, 19.0, +0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,      real(ff->f_perp(0.05, 16.0, -0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.001215, imag(ff->f_perp(0.05, 16.0, -0.5)), eps);
        }
} b_to_pi_pi_fvdv2018_form_factors_test;
