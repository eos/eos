/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#ifndef EOS_GUARD_EOS_MATHS_INTERPOLATION_HH
#define EOS_GUARD_EOS_MATHS_INTERPOLATION_HH 1

#include <vector>
#include <functional>
#include <memory>
#include <gsl/gsl_spline.h>

namespace eos
{
    class CSplineInterpolation
    {
        private:
            const std::vector<double> _data_x;
            const std::vector<double> _data_y;
            const std::function<double(const double &)> _map_x;
            const std::shared_ptr<gsl_interp_accel> _acc;
            const std::shared_ptr<gsl_interp> _interp;

        public:
            CSplineInterpolation() = delete;
            /*!
             * Stores the interpolation data and initializes gsl variables.
             *
             * @param data_x The supporting points of the x domain of the function.
             * @param data_y The corresponding function values.
             */
            CSplineInterpolation(const std::vector<double> & data_x, const std::vector<double> & data_y);

            /*!
             * Evaluate the interpolating function.
             *
             * @param x The point at which the function shall be evaluated.
             */
            double operator()(const double & x) const;
    };
}

#endif
