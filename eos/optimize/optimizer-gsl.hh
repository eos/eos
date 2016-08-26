/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_OPTIMIZE_OPTIMIZER_MINUIT2_HH
#define EOS_GUARD_EOS_OPTIMIZE_OPTIMIZER_MINUIT2_HH 1

#include <eos/optimize/optimizer.hh>

#include <gsl/gsl_multimin.h>

namespace eos
{
    class OptimizerGSL :
        public Optimizer
    {
        private:
            // target function for the minimization process
            DensityPtr _density;

            // maximum number of iterations performed to find optimium
            const unsigned _max_iterations;

            // target size for the simplex
            const double _target_size;

            // number of parameters
            const unsigned _number_of_parameters;

            // GSL vector of parameter values
            gsl_vector * _gsl_parameters;

            // GSL vector of step size per parameter;
            gsl_vector * _gsl_step_size;

            // GSL function minimization algorithm.
            const gsl_multimin_fminimizer_type * _gsl_type;

            // GSL minimization state.
            gsl_multimin_fminimizer * _gsl_state;

            // GSL function.
            gsl_multimin_function _gsl_func;

            // copy parameter values from the density to the GSL vector
            void _update_gsl_parameters(gsl_vector * gsl_parameters);

            // copy parameter values from the GSL vector to the density
            void _update_density(const gsl_vector * gsl_parameters);

            // evaluate the target function for a given GSL vector of parameters
            double _evaluate(const gsl_vector * gsl_parameters);

            // optimize until either limit of iterations has been reached, or
            // target size of the simplex has been achieved
            double _optimize();

        public:
            static double _evaluate_original_adapter(const gsl_vector * gsl_parameters, void * _this);
            static double _evaluate_negative_adapter(const gsl_vector * gsl_parameters, void * _this);

            OptimizerGSL(const DensityPtr & density, const unsigned & max_iterations, const double & target_size);

            ~OptimizerGSL();

            virtual double maximize();

            virtual double minimize();
    };
}

#endif
