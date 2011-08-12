/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#include <eos/utils/chi-squared.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <map>
#include <string>

namespace eos
{
    double
    ChiSquared::with_theory_offset(const double & theory_min, const double &, const double & theory_max,
            const double & experiment_min, const double & experiment_central, const double & experiment_max)
    {
        double sigma = std::abs(experiment_max - experiment_min);
        double chi = 0.0;

        if ((experiment_central - theory_max) > 0)
        {
            chi = (experiment_central - theory_max) / sigma;
        }
        else if ((theory_min - experiment_central) > 0)
        {
            chi = (theory_min - experiment_central) / sigma;
        }

        return chi * chi;
    }

    double
    ChiSquared::with_combined_uncertainties(const double & theory_min, const double & theory_central, const double & theory_max,
            const double & experiment_min, const double & experiment_central, const double & experiment_max)
    {
        double sigma = std::abs(experiment_max - experiment_min);
        double delta_plus = std::abs(theory_max - theory_central), delta_minus = std::abs(theory_min - theory_central);
        double chi = experiment_central - theory_central;

        if (chi > 0)
        {
            chi /= std::sqrt(power_of<2>(sigma) + power_of<2>(delta_minus));
        }
        else if (chi < 0)
        {
            chi /= std::sqrt(power_of<2>(sigma) + power_of<2>(delta_plus));
        }

        return chi * chi;
    }

    ChiSquared::Function
    ChiSquared::make(const std::string & name)
    {
        static const std::map<std::string, ChiSquared::Function> functions
        {
            std::make_pair("with-theory-offset",                &ChiSquared::with_theory_offset),
            std::make_pair("with-combined-uncertainties",       &ChiSquared::with_combined_uncertainties),
        };

        auto i = functions.find(name);
        if (functions.cend() == i)
            throw NoSuchChiSquaredError(name);

        return i->second;
    }

    NoSuchChiSquaredError::NoSuchChiSquaredError(const std::string & name) :
        Exception("No such ChiSquared function: '" + name + "'")
    {
    }
}
