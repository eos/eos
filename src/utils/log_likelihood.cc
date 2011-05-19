/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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

#include <src/utils/log_likelihood.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include <vector>
#include <algorithm>

namespace eos
{
    template <>
    struct Implementation<LogLikelihood>
    {
        // <ObservablePtr, current result of evaluation>
        std::vector<std::tuple<ObservablePtr, double>> predictions;

        // Same order as in predictions
        // <last result of evaluation>
        std::vector<double> predictions_dirty;

        // <index to predictions, min, central, max>
        std::vector<std::tuple<unsigned, double, double, double>> observations;

        // <observable name, index>
        std::map<std::string, unsigned> observable_names;

        // store which observables(index list for observations) are needed for a parameter id
        std::map<unsigned, std::vector<unsigned>> index_lists;

        Parameters parameters;

        Implementation(const Parameters & parameters) :
            parameters(parameters)
        {
        }

        static bool identical_observables(const ObservablePtr & lhs, const ObservablePtr & rhs)
        {
            // compare name
            if( lhs->name() != rhs->name())
                return false;

            // compare kinematics
            if( lhs->kinematics() != rhs->kinematics())
                return false;

            // compare options
            if (lhs->options() != rhs->options())
                return false;

            return true;
        }

        void add(const ObservablePtr & observable, const double & min, const double & central, const double & max)
        {
            if (observable->parameters() != parameters)
                throw InternalError("Likelihood::add(): Encountered observable whose parameters doesn't fit ours!");

            if ( ! index_lists.empty() )
                throw InternalError("Likelihood::add(): Add all observables before evaluating the likelihood for the first time!");

            // check if observable of that name exists already
            auto result = observable_names.insert(std::make_pair(observable->name(), predictions.size()));

            // existing observable
            if ((result.second == false) && identical_observables(observable, std::get<0>(predictions[result.first->second])))
            {
                unsigned index = result.first->second;
                observations.push_back(std::make_tuple(index, min, central, max) );
            }
            else //new observable
            {
                unsigned index = predictions.size();

                predictions.push_back(std::make_tuple(observable, std::numeric_limits<double>::quiet_NaN()));
                predictions_dirty.push_back(std::numeric_limits<double>::quiet_NaN());
                observations.push_back(std::make_tuple(index, min, central, max) );
            }
        }

        /// calculate predictions for all observables anew
        void evaluate_observables()
        {
            for (auto i = predictions.begin(), i_end = predictions.end() ; i != i_end ; ++i)
            {
                std::get<1>(*i) = std::get<0>(*i)->evaluate();
            }
        }

        void evaluate_observables(const Parameter::Id & id)
        {
            // get list of observables for this id
            auto index_list = index_lists.find(id);

            // if parameter changes for the first time, determine the  list of observables for this id
            if (index_list == index_lists.end())
            {
                std::vector<unsigned> indices;

                //loop over all observables
                unsigned index = 0;
                for (auto i = predictions.begin(),  i_end = predictions.end() ; i != i_end ; ++i, ++index)
                {
                    // search through all parameters that this observable uses
                    auto result = std::find(std::get<0>(*i)->begin(), std::get<0>(*i)->end(), id);

                    // no dependence on parameter
                    if (result == std::get<0>(*i)->end())
                       continue;

                    indices.push_back(index);
                }

                // add the index list
                index_list = index_lists.insert(std::make_pair(id, indices)).first;
            }

            // evaluate observables one by one
            for (auto i = index_list->second.begin(), i_end = index_list->second.end(); i != i_end; ++i)
            {
                // store last value of observable
                predictions_dirty.at(*i) = std::get<1>(predictions.at(*i));

                // the new value
                std::get<1>(predictions.at(*i)) = std::get<0>(predictions.at(*i))->evaluate();
            }
        }

        double operator() () const
        {
           /*
            * proper normalization
            * Could calculate only once upon initialization, but if
            * uncertainty is parameter dependent, that's ruled out.
            */
           double norm = 1.0;

           double chi_squared = 0.0;
           double sigma = 1.0;

           // compare predictions with observations
           for (auto i = observations.cbegin(), i_end = observations.cend() ; i != i_end ; ++i)
           {
               // get prediction
               double value = std::get<1>(predictions[std::get<0>(*i)]);

               // get most likely observed value
               double central = std::get<2>(*i);

               // allow for asymmetric Gaussian uncertainty
               if (value > central)
                   sigma = std::get<3>(*i) - central;
               else
                   sigma = central - std::get<1>(*i);

               double chi = (value - central) / sigma;

               chi_squared += chi * chi;

               norm /= std::sqrt(2.0 * M_PI) * sigma;
           }

           return std::log(norm) - chi_squared / 2.0;
        }
    };

    LogLikelihood::LogLikelihood(const Parameters & parameters) :
        PrivateImplementationPattern<LogLikelihood>(new Implementation<LogLikelihood>(parameters))
    {
    }

    LogLikelihood::~LogLikelihood()
    {
    }

    void
    LogLikelihood::add(const ObservablePtr & observable, const double & min, const double & central, const double & max)
    {
        _imp->add(observable, min, central, max);
    }

    LogLikelihood
    LogLikelihood::clone() const
    {
        Parameters parameters(_imp->parameters.clone());
        LogLikelihood result(parameters);

        for (auto o = _imp->observations.cbegin(), o_end = _imp->observations.cend() ; o != o_end ; ++o)
        {
            // create clone with independent parameters object
            ObservablePtr ptr = std::get<0>(_imp->predictions.at(std::get<0>(*o)))->clone(parameters);
            result.add(ptr, std::get<1>(*o), std::get<2>(*o), std::get<3>(*o));
        }

        return result;
    }

    Parameters
    LogLikelihood::parameters() const
    {
        return _imp->parameters;
    }

    const std::vector<std::tuple<ObservablePtr, double>> &
    LogLikelihood::predictions() const
    {
        return _imp->predictions;
    }

    double
    LogLikelihood::operator() ()
    {
        _imp->evaluate_observables();
        return (*_imp)();
    }

    double
    LogLikelihood::operator() (const Parameter::Id& id)
    {
        _imp->evaluate_observables(id);
        return (*_imp)();
    }

    void
    LogLikelihood::reset()
    {
        for (unsigned i = 0 ; i < _imp->predictions.size() ; ++i)
        {
            std::get<1>(_imp->predictions[i]) = _imp->predictions_dirty[i];
        }
    }
}
