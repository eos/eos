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

#include <eos/optimize/optimizer-gsl.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    OptimizerGSL::OptimizerGSL(const DensityPtr & density, const unsigned & max_iterations, const double & target_size) :
        _density(density),
        _max_iterations(max_iterations),
        _target_size(target_size),
        _number_of_parameters(std::distance(density->begin(), density->end())),
        _gsl_parameters(gsl_vector_alloc(_number_of_parameters)),
        _gsl_step_size(gsl_vector_alloc(_number_of_parameters)),
        _gsl_type(gsl_multimin_fminimizer_nmsimplex2),
        _gsl_state(gsl_multimin_fminimizer_alloc(_gsl_type, _number_of_parameters))
    {
        _update_gsl_parameters(_gsl_parameters);

        unsigned i = 0;
        for (auto p = density->begin(), p_end = density->end() ; p != p_end ; ++p, ++i)
        {
            gsl_vector_set(_gsl_step_size, i, (p->max - p->min) / 100);
        }
    }

    OptimizerGSL::~OptimizerGSL()
    {
        gsl_vector_free(_gsl_step_size);
        gsl_vector_free(_gsl_parameters);
        gsl_multimin_fminimizer_free(_gsl_state);
    }

    void
    OptimizerGSL::_update_gsl_parameters(gsl_vector * gsl_parameters)
    {
        unsigned i = 0;

        for (auto p = _density->begin(), p_end = _density->end() ; p != p_end ; ++p, ++i)
        {
            gsl_vector_set(gsl_parameters, i, p->parameter->evaluate());
        }
    }

    void
    OptimizerGSL::_update_density(const gsl_vector * gsl_parameters)
    {
        unsigned i = 0;

        for (auto p = _density->begin(), p_end = _density->end() ; p != p_end ; ++p, ++i)
        {
            p->parameter->set(gsl_vector_get(gsl_parameters, i));
        }
    }

    double
    OptimizerGSL::_evaluate(const gsl_vector * gsl_parameters)
    {
        _update_density(gsl_parameters);

        return _density->evaluate();
    }

    double
    OptimizerGSL::_evaluate_original_adapter(const gsl_vector * gsl_parameters, void * _this)
    {
        return static_cast<OptimizerGSL *>(_this)->_evaluate(gsl_parameters);
    }

    double
    OptimizerGSL::_evaluate_negative_adapter(const gsl_vector * gsl_parameters, void * _this)
    {
        return -static_cast<OptimizerGSL *>(_this)->_evaluate(gsl_parameters);
    }

    double
    OptimizerGSL::_optimize()
    {
        unsigned iterations = 0;
        int status;

        do
        {
            iterations++;
            status = gsl_multimin_fminimizer_iterate(_gsl_state);

            if (status)
                break;

            double simplex_size = gsl_multimin_fminimizer_size(_gsl_state);
            status = gsl_multimin_test_size(simplex_size, _target_size);

            if (GSL_SUCCESS == status)
            {
                return gsl_multimin_fminimizer_minimum(_gsl_state);
            }

        }
        while ((GSL_CONTINUE == status) && (iterations < _max_iterations));

        if (GSL_SUCCESS == status)
        {
            return gsl_multimin_fminimizer_minimum(_gsl_state);
        }

        throw OptimizerError("GSL multimin did not converge after " + stringify(_max_iterations) + " iterations!");
    }

    double
    OptimizerGSL::maximize()
    {
        _gsl_func.n = _number_of_parameters;
        _gsl_func.f = &OptimizerGSL::_evaluate_negative_adapter;
        _gsl_func.params = static_cast<void *>(this);

        gsl_multimin_fminimizer_set(_gsl_state, &_gsl_func, _gsl_parameters, _gsl_step_size);

        return _optimize();
    }

    double
    OptimizerGSL::minimize()
    {
        _gsl_func.n = _number_of_parameters;
        _gsl_func.f = &OptimizerGSL::_evaluate_original_adapter;
        _gsl_func.params = static_cast<void *>(this);

        gsl_multimin_fminimizer_set(_gsl_state, &_gsl_func, _gsl_parameters, _gsl_step_size);

        return _optimize();
    }
}
