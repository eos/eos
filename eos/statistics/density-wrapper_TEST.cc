/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Frederik Beaujean
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

#include <eos/statistics/density-wrapper_TEST.hh>
#include <eos/statistics/simple-parameters.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>
#include <test/test.hh>

#include <cmath>
#include <map>

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
}

class DensityTest :
    public TestCase
{
    public:
        DensityTest() :
            TestCase("density_wrapper_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-13;

            // create
            {
                DensityWrapper::WrappedDensity wrapped_density(::multivariate_unit_normal_pdf);
                DensityWrapper density(wrapped_density);

                // access via name
                density.add_parameter("x", -5, 5);
                density.parameters()["x"] = 1.5;

                // access via iterator
                density.add_parameter("y", -5, 5);
                *(++density.begin())->parameter = -0.3;

                static const double result = -3.0078770664093453;
                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), result, eps);

                // copy
                DensityWrapper density_copy = density;
                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), result, eps);

                // clone
                DensityPtr density_clone = density.clone();
                TEST_CHECK_RELATIVE_ERROR(density_clone->evaluate(), result, eps);
            }

            // modify
            {
                DensityWrapper density(make_multivariate_unit_normal(2));
                density.begin()->parameter->set(1.5);
                (++density.begin())->parameter->set(-0.3);

                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), -3.0078770664093453, eps);
            }
        }
} density_wrapper_test;
