/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016 Danny van Dyk
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

#include <eos/optimize/optimizer.hh>
#include <eos/optimize/optimizer-gsl.hh>
#include <eos/statistics/density-wrapper.hh>
#include <eos/statistics/simple-parameters.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>
#include <test/test.hh>

#include <cmath>
#include <map>

#include <iostream>

using namespace test;
using namespace eos;

namespace
{
    /*!
     * A multivariate normal distribution with
     * zero mean and unit covariance on the log scale
     */
    double multivariate_unit_normal_pdf(const std::vector<double> parameters)
    {
        double result = 0;

        for (const auto & p : parameters)
        {
            result += -log(std::sqrt(2 * M_PI)) - power_of<2>(p) / 2.0;
        }

        return result;
    }

    /*!
     * A multivariate normal distribution with
     * mean mu_i = i mod 5, and unit covariance on the log scale
     */
    double shifted_multivariate_unit_normal_pdf(const std::vector<double> parameters)
    {
        double result = 0;
        unsigned i = 0;

        for (const auto & p : parameters)
        {
            double mean = (i % 5) * 1.0 + 0.0;
            result += /*-log(std::sqrt(2 * M_PI))*/ - power_of<2>(p - mean) / 2.0;

            ++i;
        }

        return result;
    }
}

namespace eos
{
     DensityWrapper
     make_multivariate_unit_normal(const unsigned & ndim)
     {
         DensityWrapper::WrappedDensity wrapped_density(::multivariate_unit_normal_pdf);
         DensityWrapper density(wrapped_density);

         for (unsigned i = 0 ; i < ndim ; ++i)
         {
             density.add_parameter(std::string("par") + stringify(i), -5, 5);
         }

         return density;
     }

     DensityWrapper
     make_shifted_multivariate_unit_normal(const unsigned & ndim)
     {
         DensityWrapper::WrappedDensity wrapped_density(::shifted_multivariate_unit_normal_pdf);
         DensityWrapper density(wrapped_density);

         for (unsigned i = 0 ; i < ndim ; ++i)
         {
             density.add_parameter(std::string("par") + stringify(i), -5, 5);
         }

         return density;
     }
}

class OptimizerGSLTest :
    public TestCase
{
    public:
        OptimizerGSLTest() :
            TestCase("optimizer_gsl_test")
        {
        }

        virtual void run() const
        {
            // multivariate_unit_normal: 3D
            {
                static const double eps = 1e-5;
                static const double target_size = 1e-6;
                static const unsigned max_iterations = 200;

                DensityPtr density(new DensityWrapper(make_multivariate_unit_normal(3)));

                OptimizerPtr optimizer(new OptimizerGSL(density, max_iterations, target_size));
                optimizer->maximize();

                for (auto & p : *density)
                {
                    TEST_CHECK_NEARLY_EQUAL(p.parameter->evaluate(), 0.0, eps);
                }
            }

            // multivariate_unit_normal: 20D
            {
                static const double eps = 1e-5;
                static const double target_size = 1e-6;
                static const unsigned max_iterations = 2000;

                DensityPtr density(new DensityWrapper(make_multivariate_unit_normal(20)));

                OptimizerPtr optimizer(new OptimizerGSL(density, max_iterations, target_size));
                optimizer->maximize();

                for (auto & p : *density)
                {
                    TEST_CHECK_NEARLY_EQUAL(p.parameter->evaluate(), 0.0, eps);
                }
            }

            // shifted_multivariate_unit_normal: 5D
            {
                static const double eps = 1e-5;
                static const double target_size = 1e-6;
                static const unsigned max_iterations = 500;

                DensityPtr density(new DensityWrapper(make_shifted_multivariate_unit_normal(5)));

                OptimizerPtr optimizer(new OptimizerGSL(density, max_iterations, target_size));
                optimizer->maximize();

                unsigned i = 0;
                for (auto & p : *density)
                {
                    unsigned index = i++;
                    double value = index;

                    TEST_CHECK_NEARLY_EQUAL(p.parameter->evaluate(), value, eps);
                }
            }

            // shifted_multivariate_unit_normal: 23D
            {
                static const double eps = 1e-5;
                static const double target_size = 1e-7;
                static const unsigned max_iterations = 10000;

                DensityPtr density(new DensityWrapper(make_shifted_multivariate_unit_normal(20)));

                unsigned i = 0;

                // for this test, the simplex optimization needs a good initial
                // point.
                for (auto & p : *density)
                {
                    unsigned index = i++;
                    double value = (index % 5) * 1.01 + (index % 2 == 0 ? -1.0 : +1.0) * 0.05;

                    p.parameter->set(value);
                }

                OptimizerPtr optimizer(new OptimizerGSL(density, max_iterations, target_size));
                optimizer->maximize();

                for (auto & p : *density)
                {
                    unsigned index = i++;
                    double value = index % 5;

                    TEST_CHECK_NEARLY_EQUAL(p.parameter->evaluate(), value, eps);
                }
            }
        }
} optimizer_gsl_test;
