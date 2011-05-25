/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <src/utils/analysis.hh>
#include <src/utils/log.hh>
#include <src/utils/markov_chain.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/stringify.hh>
#include <cmath>
#include <limits>
#include <numeric>
#include <map>

#include <gsl/gsl_rng.h>

namespace eos
{
    template<>
    struct Implementation<MarkovChain>
    {
        // was the last proposed move accepted?
        bool accept_proposal;

        // The Analysis that we are conducting.
        std::shared_ptr<Analysis> analysis;

        // all the info about parameters that is needed to setup sampling
        std::vector<ParameterDescription> parameter_descriptions;

        // get sample values directly from the prior for discrete parameters
        std::map<unsigned, LogPriorPtr> discrete_priors;

        // Random number generator, unique to this chain
        gsl_rng * rng;

        // scale factors used in proposal functions
        std::vector<double> scales;

        // how far are we in the current run?
        // Cleared after each reset, e.g. in prerun
        unsigned current_iteration;

        // index of current parameter that is changed
        unsigned current_parameter;

        // info for our current point
        MarkovChain::InfoAtPoint current;

        // history of the random walk
        MarkovChain::History history;

        // proposed value of current parameter
        double prop;

        // info for our proposed point
        MarkovChain::InfoAtPoint proposal;

        // total number of iterations in this/the last sampling run
        unsigned run_iterations;

        // overall statistics
        MarkovChain::Stats stats;

        // sample variance of param values (Welford's method)
        std::vector<double> variance_of_par_temp;

        // sample variance of log(posterior) (Welford's method)
        double variance_of_posterior_temp;

        Implementation(const std::shared_ptr<Analysis> & analysis, unsigned long seed) :
            analysis(analysis)
        {
            // setup Mersenne-Twister RN generator using default seed
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, seed);

            initialize();
        }

        ~Implementation<MarkovChain>()
        {
            // free RN generator
            gsl_rng_free(rng);
        }

        // clear this chains' history
        void clear()
        {
            history.points.clear();
            history.log_likelihood.clear();
            history.log_posterior.clear();
            history.log_prior.clear();
        }

        // dump history of this chain into a ScanFile's DataSet
        void dump_history(ScanFile::DataSet & data_set)
        {
            Log::instance()->message("markov_chain.dump_history", ll_debug)
                << "Dumping " << history.points.size() << " records";

            ScanFile::WriteBuffer buffer(data_set.fields(), history.points.size());

            std::vector<double> tuple;

            for (unsigned i = 0 ; i < history.points.size() ; ++i)
            {
                tuple = history.points[i];
                tuple.push_back(history.log_posterior[i]);

                buffer << tuple;

                // once the buffer is full, write to disk
                if (buffer.capacity() == buffer.size())
                {
                    data_set << buffer;
                    buffer.clear();
                }
            }

            // write the remainder to disk
            data_set << buffer;
            buffer.clear();
        }

        /*!
         * calculate posterior etc at the proposal point
         */
        void evaluate_point()
        {
            // check proposal point
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                if ((proposal.point[i] < parameter_descriptions[i].min)
                        || (proposal.point[i] > parameter_descriptions[i].max))
                {
                    throw InternalError("MarkovChain::evaluate_point: parameter '" + parameter_descriptions[i].parameter.name()
                            + "' = " + stringify(proposal.point[i]) + " not in valid range ["
                            + stringify(parameter_descriptions[i].min) + "," + stringify(parameter_descriptions[i].max) + "]"
                            + " in iteration " + stringify(current_iteration));
                }
            }

            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                if (parameter_descriptions[i].parameter != current.point[i] )
                    throw InternalError("MarkovChain::evaluate_point: parameter '" + parameter_descriptions[i].parameter.name()
                            + "' = " + stringify(parameter_descriptions[i].parameter()) + " doesn't match current point!");
            }

            // change Parameter object
            parameter_descriptions[current_parameter].parameter = prop;

            // let likelihood only evaluate observables that use the parameter we just changed
            proposal.log_likelihood = (*analysis->log_likelihood())(parameter_descriptions[current_parameter].parameter.id());
            proposal.log_prior = analysis->log_prior();
            proposal.log_posterior = proposal.log_prior + proposal.log_likelihood;
        }

        // called from ctor only at beginning
        void initialize()
        {
            // copy the information about parameters, their ranges and
            // whether they are nuisance parameters or not
            parameter_descriptions = analysis->parameter_descriptions();

            //find discrete parameters
            unsigned index = 0;
            for (auto i = parameter_descriptions.begin(), i_end = parameter_descriptions.end(); i != i_end; ++i, ++index)
            {
                if (i->discrete)
                    discrete_priors.insert(std::make_pair(index, analysis->log_prior(i->parameter.name())));
            }

            // initialize statistics
            reset(true);

            // now update storage capacities
            scales.resize(parameter_descriptions.size(), 1);
            current.point.resize(parameter_descriptions.size(), 0);
            proposal.point.resize(parameter_descriptions.size(), 0);

            // by default save points and posterior values
            history.keep = true;

            // set initial points
            initial_points();
        }

        /*
         * Choose random starting positions for each chain
         */
        void initial_points()
        {
            // random starting point, equal current and proposal.
            //   x_{init} = x_{min} + U * (x_{max}-x_{min})
            unsigned i = 0;
            for (auto p = parameter_descriptions.begin(), p_end = parameter_descriptions.end(); p != p_end; ++p,++i)
            {
                if(p->discrete)
                    current.point[i] = discrete_priors[i]->sample(rng);
                else
                    current.point[i] = p->min + uniform_random_number() * ( p->max - p->min);

                // update Parameters object
                proposal.point[i] = current.point[i];
                p->parameter = current.point[i];
            }

            // reuse code. Pretend we propose to change first param. to its initial value
            if (parameter_descriptions.size() > 0)
            {
                current_parameter = 0;
                prop = current.point[current_parameter];

                // evaluate likelihood without argument, so all observables are calculated once
                proposal.log_likelihood = (*analysis->log_likelihood())();
                proposal.log_prior = analysis->log_prior();
                proposal.log_posterior = proposal.log_prior + proposal.log_likelihood;
                move();

                // don't update!
            }

            // setup statistics
            stats.mode_of_posterior = current.log_posterior;
            stats.parameters_at_mode = current.point;
        }

        /*
         * Propose point, set it up and check if the point is in range.
         *
         * @return true, if proposal point is in range.
         */
        bool new_proposal_point()
        {
            ParameterDescription& descr = parameter_descriptions[current_parameter];

            // ensure that we start from current point
            // undo any changes for previous parameter.
            // to cover the case where parameter = 0, need modulus to wrap around
            unsigned previous_par_index = (current_parameter == 0) ? parameter_descriptions.size() - 1 : current_parameter - 1;
            proposal.point[previous_par_index] = current.point[previous_par_index];

            if (descr.discrete)
                prop = discrete_priors[current_parameter]->sample(rng);
            else
            {
                double scale = scales[current_parameter];

                // u in [0,1[
                double u = uniform_random_number();

                // inverse transform method for Cauchy variable around zero
                double x = scale * std::tan(M_PI * (u - 1 / 2.)) + 0;

                // proposal = current point plus move
                double par_range = descr.max - descr.min;
                prop = current.point[current_parameter] + x * par_range;
            }

            // save point and update parameters
            proposal.point[current_parameter] = prop;

            // check if point outside valid range. Then it is rejected and no likelihood
            // evaluation is needed
            if ((prop < descr.min) || (prop > descr.max))
            {
                return false;
            }

            return true;
        }

        // copy proposal value once move is accepted */
        inline void move()
        {
            current.log_likelihood = proposal.log_likelihood;
            current.log_prior = proposal.log_prior;
            current.log_posterior = proposal.log_posterior;
            current.point[current_parameter] = prop;
        }

        /** Clear all statistics and counter.
         *
         * Do not change current position nor scale.
         * @par hard loose all statistics (use after prerun before main run)
         * */
        void reset(bool hard = false)
        {
            current_iteration = 0;

            // clear statistics
            stats.iterations_accepted.clear();
            stats.iterations_accepted.resize(parameter_descriptions.size(), 0.0);
            stats.iterations_rejected.clear();
            stats.iterations_rejected.resize(parameter_descriptions.size(), 0.0);

            if (hard)
            {
                stats.iterations_total = 0;

                stats.mean_of_parameters.clear();
                stats.mean_of_parameters.resize(parameter_descriptions.size(), 0.0);
                stats.mean_of_posterior = 0.0;

                stats.variance_of_parameters.clear();
                stats.variance_of_parameters.resize(parameter_descriptions.size(), 0.0);
                variance_of_par_temp.clear();
                variance_of_par_temp.resize(parameter_descriptions.size(), 0.0);
                stats.variance_of_posterior = 0.0;
                variance_of_posterior_temp = 0.0;

                stats.mode_of_posterior = -std::numeric_limits<double>::max();
            }
        }

        // undo changes to Parameter object
        inline void revert()
        {
            parameter_descriptions[current_parameter].parameter = current.point[current_parameter];
            // reload old observable values
            analysis->log_likelihood()->reset();
        }

        // set the number of iterations for next run and go
        void run(unsigned iterations)
        {
            Log::instance()->message("markov_chain.run", ll_debug)
                << "Running " << iterations << " iterations";

            reset();

            // make sure everything is fine __before__ we start
            self_check();

            // loop over iterations
            for (current_iteration = 0 ; current_iteration < iterations ; ++current_iteration)
            {
                // loop over parameters, one at a time
                for (current_parameter = 0 ; current_parameter < parameter_descriptions.size() ; ++current_parameter)
                {
                    accept_proposal = false;

                    // proposal point updated!
                    bool in_range = new_proposal_point();

                    if (in_range)
                    {
                        // calc. likelihood, log_prior. Parameter object updated!
                        evaluate_point();

                        // compare, throw dice
                        double r = std::log(uniform_random_number());

                        if (r < proposal.log_posterior - current.log_posterior)
                        {
                            accept_proposal = true;
                            move();
                        }
                        else
                        {
                            // restore previous state of Parameters
                            revert();
                        }
                    }

                    // save points, update statistics etc
                    update();
                }

            }

            // we are done. store how many iterations we had in total
            stats.iterations_total += iterations;
            run_iterations = iterations;
        }

        // check consistency of configuration, throw exception
        void self_check()
        {
            if (parameter_descriptions.size() == 0)
            {
                throw InternalError("MarkovChain::selfCheck(): Number of parameters does not exceed 0");
            }
        }

        // save points, update statistics
        void update()
        {
            // store points
            if (history.keep)
            {
                history.points.push_back(current.point);
                history.log_likelihood.push_back(current.log_likelihood);
                history.log_posterior.push_back(current.log_posterior);
                history.log_prior.push_back(current.log_prior);
            }

            if (accept_proposal)
            {
                ++stats.iterations_accepted[current_parameter];
            }
            else
            {
                ++stats.iterations_rejected[current_parameter];
            }

            // count iterations for this parameter since (pre|main) run started. start index at 0, so need +1
            double total_iterations_since_reset = stats.iterations_total + (current_iteration + 1.0);

            if (current.log_posterior > stats.mode_of_posterior)
            {
                stats.mode_of_posterior = current.log_posterior;
                stats.parameters_at_mode = current.point;
            }

            // update mean, keep copy for variance calculation below
            double tmp_parameter = stats.mean_of_parameters[current_parameter];
            stats.mean_of_parameters[current_parameter] += (current.point[current_parameter] - tmp_parameter) / total_iterations_since_reset;
            double tmpPost = stats.mean_of_posterior;
            stats.mean_of_posterior += (current.log_posterior - tmpPost) / total_iterations_since_reset;

            // update variance using Welford's method
            // see http://www.johndcook.com/standard_deviation.html,
            // Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd edition
            if (total_iterations_since_reset < 2)
            {
                variance_of_par_temp[current_parameter] = 0;
                variance_of_posterior_temp = 0;
            }
            else
            {
                variance_of_par_temp[current_parameter] += (current.point[current_parameter] - tmp_parameter) *
                    (current.point[current_parameter] - stats.mean_of_parameters[current_parameter]);
                variance_of_posterior_temp += (current.log_posterior - tmpPost) * (current.log_posterior - stats.mean_of_posterior);

                stats.variance_of_parameters[current_parameter] = variance_of_par_temp[current_parameter] / (total_iterations_since_reset - 1);
                stats.variance_of_posterior = variance_of_posterior_temp / (total_iterations_since_reset - 1);
            }
        }

        // return random number in [0,1[
        inline double uniform_random_number() { return gsl_rng_uniform(rng); }
    };

    MarkovChain::MarkovChain(const std::shared_ptr<Analysis> & analysis, unsigned long seed) :
        PrivateImplementationPattern<MarkovChain>(new Implementation<MarkovChain>(analysis, seed))
    {
    }

    MarkovChain::~MarkovChain()
    {
    }

    void
    MarkovChain::clear()
    {
        _imp->clear();
    }

    void
    MarkovChain::dump_history(ScanFile::DataSet & data_set)
    {
        _imp->dump_history(data_set);
    }

    double
    MarkovChain::get_scale(const unsigned & index) const
    {
        return _imp->scales[index];
    }

    void
    MarkovChain::set_scale(const double & scale)
    {
        for (auto s = _imp->scales.begin(), s_end = _imp->scales.end() ; s != s_end ; ++s)
        {
            *s = scale;
        }
    }

    void
    MarkovChain::set_scale(const unsigned & index, const double & scale)
    {
        _imp->scales[index] = scale;
    }

    const MarkovChain::InfoAtPoint &
    MarkovChain::info_at_current() const
    {
        return _imp->current;
    }

    const MarkovChain::InfoAtPoint &
    MarkovChain::info_at_proposal() const
    {
        return _imp->proposal;
    }

    const unsigned &
    MarkovChain::iterations_last_run() const
    {
        return _imp->run_iterations;
    }

    void
    MarkovChain::keep_history(bool keep)
    {
        _imp->history.keep = keep;
    }

    const std::vector<ParameterDescription> &
    MarkovChain::parameter_descriptions() const
    {
        return _imp->parameter_descriptions;
    }

    bool
    MarkovChain::proposal_accepted() const
    {
        return _imp->accept_proposal;
    }

    void
    MarkovChain::reset(bool hard)
    {
        _imp->reset(hard);
    }

    void
    MarkovChain::run(const unsigned & iterations)
    {
        _imp->run(iterations);
    }

    const MarkovChain::Stats &
    MarkovChain::statistics() const
    {
        return _imp->stats;
    }

    std::ostream &
    operator<< (std::ostream & lhs, const MarkovChain::InfoAtPoint & rhs)
    {
        lhs << "point = ( ";

        for (auto p = rhs.point.cbegin(), p_end = rhs.point.cend() ; p != p_end ; ++p)
        {
            lhs << *p << " ";
        }

        lhs << "), prior = " << rhs.log_prior << ", likelihood = " << rhs.log_likelihood << ", posterior = " << rhs.log_posterior;

        return lhs;
    }
}
