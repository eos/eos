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

#include <interpolation.hh>
#include <gsl/gsl_errno.h>
#include <eos/utils/exception.hh>

namespace eos
{
    CSplineInterpolation::CSplineInterpolation(const std::vector<double> & data_x, const std::vector<double> & data_y):
        _data_x(data_x),
        _data_y(data_y),
        _acc(gsl_interp_accel_alloc(), &gsl_interp_accel_free),
        _interp(gsl_interp_alloc(gsl_interp_cspline, _data_x.size()), &gsl_interp_free)
    {
        if (_data_x.size() != _data_y.size())
        {
            throw InternalError("Interpolation: dimensions of x and y data does not match");
        }
        gsl_interp_init(_interp.get(), &_data_x[0], &_data_y[0], _data_x.size());
    }

    double CSplineInterpolation::operator()(const double & x) const
    {
        double res = 0;
        int gsl_status = 0;

        gsl_status = gsl_interp_eval_e(_interp.get(), &_data_x[0], &_data_y[0], x, _acc.get(), &res);
        if (gsl_status)
        {
            throw eos::GSLError(gsl_strerror(gsl_status));
        }

        return res;
    }
}
