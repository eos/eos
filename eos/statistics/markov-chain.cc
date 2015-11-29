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

#include <config.h>

#include <eos/statistics/markov-chain.hh>
#include <eos/statistics/proposal-functions.hh>
#include <eos/utils/density.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>

#include <gsl/gsl_rng.h>

namespace eos
{
    template<>
    struct Implementation<MarkovChain>
    {
        // The Analysis that we are conducting.
        DensityPtr density;

        // The proposal function
        std::shared_ptr<MarkovChain::ProposalFunction> proposal_function;

        // all the info about parameters that is needed to setup sampling
        std::vector<ParameterDescription> parameter_descriptions;

        // Random number generator, unique to this chain
        gsl_rng * rng;

        // was the last proposed move accepted?
        bool accept_proposal;

        // how far are we in the current run?
        // Cleared after each reset, e.g. in prerun
        unsigned current_iteration;

        // info for our current point
        MarkovChain::State current;

        // info for our proposed point
        MarkovChain::State proposal;

        // history of the random walk
        MarkovChain::History history;

        // total number of iterations in this/the last sampling run
        unsigned run_iterations;

        // overall statistics
        MarkovChain::Stats stats;

        // sample variance of param values (Welford's method)
        std::vector<double> welford_data_parameters;

        // sample variance of log(density) (Welford's method)
        double welford_data_density;

        // Output data types
        typedef hdf5::Array<1, double> SampleType;
        const SampleType sample_type;

        Implementation(const DensityPtr & density, unsigned long seed, const std::shared_ptr<MarkovChain::ProposalFunction> & proposal_function) :
            density(density->clone()),
            sample_type
            {
                "samples",
                { std::distance(density->begin(), density->end()) + 1ul },
            }

        {
            if (! proposal_function)
                throw InternalError("MarkovChain needs a non-empty proposal function");
            this->proposal_function = proposal_function->clone(),

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
            history.states.clear();
        }

        void dump_history(hdf5::File & file, const std::string & data_set_base_name, const unsigned & last_iterations) const
        {

            /* store samples */

            if (history.states.size() < last_iterations)
                throw InternalError("MarkovChain::dump_history: Cannot store more samples (" + stringify(last_iterations)
                    + ") than there are in history (" + stringify(history.states.size()) + ").");

            unsigned sample_record_length = parameter_descriptions.size() + 1;

            // we could get into trouble if we attempt to create a data set a 2nd time
            auto data_set = file.create_or_open_data_set(data_set_base_name + "/samples", sample_type);

            // parameter values + density
            std::vector<double> record(sample_record_length);
            for (auto s = history.states.cend() - last_iterations, s_end = history.states.cend() ; s != s_end ; ++s)
            {
                std::copy(s->point.cbegin(), s->point.cend(), record.begin());
                record.back() = s->log_density;
                data_set << record;
            }

            /* store (mode, max log(density) */

            auto data_set_mode = file.create_or_open_data_set(data_set_base_name + "/stats/mode", sample_type);
            std::copy(stats.parameters_at_mode.cbegin(), stats.parameters_at_mode.cend(), record.begin());
            record.back() = stats.mode;
            data_set_mode << record;
        }

        // store proposal density state
        void dump_proposal(hdf5::File & file, const std::string & data_set_base_name) const
        {
            proposal_function->dump_state(file, data_set_base_name + "/proposal");
        }

        // calculate density etc at the proposal point
        void evaluate_proposal()
        {
            //todo this is for debug purposes, and should never throw during production run
#if 1
            // check proposal point
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                if ((proposal.point[i] < parameter_descriptions[i].min) ||
                    (proposal.point[i] > parameter_descriptions[i].max))
                {
                    throw InternalError("MarkovChain::evaluate_point: parameter '" + parameter_descriptions[i].parameter->name()
                            + "' = " + stringify(proposal.point[i]) + " not in valid range ["
                            + stringify(parameter_descriptions[i].min) + "," + stringify(parameter_descriptions[i].max) + "]"
                            + " in iteration " + stringify(current_iteration));
                }
            }

            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                if (parameter_descriptions[i].parameter->evaluate() != current.point[i] )
                    throw InternalError("MarkovChain::evaluate_point: parameter '" + parameter_descriptions[i].parameter->name()
                            + "' = " + stringify(parameter_descriptions[i].parameter->evaluate())
                            + " doesn't match current point " + stringify(current.point[i])
                            + " in iteration " + stringify(current_iteration)
                            + ". Check if thread safety is violated due to incorrect ParameterDescription cloning");
            }
#endif
            // change Parameter object
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                parameter_descriptions[i].parameter->set(proposal.point[i]);
            }

            // finally evaluate the target density
            proposal.log_density = density->evaluate();
        }

        // called from ctor only at beginning
        void initialize()
        {
            // copy the information about parameters, their ranges and
            // whether they are nuisance parameters or not
            parameter_descriptions.clear();
            std::copy(density->begin(), density->end(), std::back_inserter(parameter_descriptions));

            // initialize statistics
            reset(true);

            // now update storage capacities
            current.point.resize(parameter_descriptions.size(), 0);
            proposal.point.resize(parameter_descriptions.size(), 0);

            // by default save points and density values
            history.keep = true;

            // uniformly distributed random starting point
            //   x_{init} = x_{min} + U * (x_{max}-x_{min})
            auto i = current.point.begin();
            for (auto p = parameter_descriptions.begin(), p_end = parameter_descriptions.end() ; p != p_end ; ++p, ++i)
            {
                //  don't draw from priors: they don't know about restricted ranges
                //TODO: draw from prior for analysis? No, but support setting initial point
                const double value = p->min + uniform_random_number() * (p->max - p->min);
                *i = value;
                p->parameter->set(value);
            }

            // need the density value at initial
            current.log_density = density->evaluate();

            // set proposal to current
            proposal = current;

            Log::instance()->message("markov_chain.ctor", ll_debug)
                << "Starting chain at: " << current;

            // setup mode
            stats.mode = current.log_density;
            stats.parameters_at_mode = current.point;
        }

        /*
         * Propose point, set it up and check if the point is in range.
         *
         * @return true, if proposal point is in range.
         */
        bool accept()
        {
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                // check if proposed point is outside valid range. Then it is rejected and no likelihood
                // evaluation is needed

                if ((proposal.point[i] < parameter_descriptions[i].min) || (proposal.point[i] > parameter_descriptions[i].max))
                {
                    stats.iterations_invalid++;
                    return false;
                }
            }

            // evaluate density at proposal point
            evaluate_proposal();

            // compute the Metropolis-Hastings factor
            double log_u = std::log(uniform_random_number());
            double log_r_post = proposal.log_density - current.log_density;
            double log_r_prop = proposal_function->evaluate(current, proposal) - proposal_function->evaluate(proposal, current);
            double log_r = log_r_post + log_r_prop;

            if ( ! std::isfinite(log_r))
                throw InternalError("MarkovChain::run: isfinite failed, either from a bad density value ("
                                    + stringify(log_r_post, 6) +
                                    ") or (more likely) from a bad value in the proposal evaluation ("
                                    + stringify(log_r_prop, 6) +
                                    "). Check if proposal covariance matrix is not invertible");

            if (log_u < log_r)
                return true;

            return false;
        }

        // copy proposal value once move is accepted */
        inline void move()
        {
            current = proposal;
        }

        static void read_history(hdf5::File & file, const std::string & data_set_base_name,
                                 const unsigned & dimension, MarkovChain::History & history)
        {
            // get meta info
            SampleType sample_type
            {
                "samples",
                { dimension + 1 },
            };
            auto data_set = file.open_data_set(data_set_base_name + "/samples", sample_type);
            std::vector<double> record(dimension + 1);
            MarkovChain::State state;
            state.point.resize(dimension);
            for (unsigned i = 0 ; i < data_set.records() ; ++i)
            {
                data_set >> record;

                std::copy(record.begin(), record.end() - 1, state.point.begin());
                state.log_density = record.back();
                history.states.push_back(state);
            }
        }

        static void read_proposal(hdf5::File & file, const std::string & data_set_base_name,
                                  const std::string & proposal_name, const unsigned & dimension,
                                  ProposalFunctionPtr & proposal)
        {
            proposal = proposal_functions::Factory::make(file, data_set_base_name, proposal_name, dimension);
        }

        static void read_stats(hdf5::File & file, const std::string & data_set_base_name,
                               const unsigned & dimension, MarkovChain::Stats & stats)
        {
            SampleType sample_type
            {
                "samples",
                { dimension + 1 },
            };
            std::vector<double> record(dimension + 1);
            auto data_set_mode = file.open_data_set(data_set_base_name + "/stats/mode", sample_type);
            // jump to last
            data_set_mode.end();
            data_set_mode >> record;

            stats.mode = record.back();
            stats.parameters_at_mode.resize(dimension);
            std::copy(record.begin(), record.end() - 1, stats.parameters_at_mode.begin());

            //todo remove later when bug is fixed
            if (record.front() == 0 && record.back() == 0)
            {
                Log::instance()->message("MarkovChain::read_stats", ll_informational)
                    << "Using next to last record for the mode, as last record seems invalid";
                data_set_mode.set_index(data_set_mode.records() - 2);
                data_set_mode >> record;
                std::copy(record.begin(), record.end() - 1, stats.parameters_at_mode.begin());
            }
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
            stats.iterations_accepted = 0;
            stats.iterations_rejected = 0;
            stats.iterations_invalid = 0;

            if (hard)
            {
                stats.iterations_total = 0;

                stats.mean_of_parameters.clear();
                stats.mean_of_parameters.resize(parameter_descriptions.size(), 0.0);
                stats.mean_of_log_density = 0.0;

                stats.variance_of_parameters.clear();
                stats.variance_of_parameters.resize(parameter_descriptions.size(), 0.0);
                welford_data_parameters.clear();
                welford_data_parameters.resize(parameter_descriptions.size(), 0.0);
                stats.variance_of_log_density = 0.0;
                welford_data_density = 0.0;

                stats.mode = -std::numeric_limits<double>::max();
            }
        }

        // undo changes to Parameter object
        inline void revert()
        {
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                parameter_descriptions[i].parameter->set(current.point[i]);
            }
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
                proposal_function->propose(proposal, current, rng);

                accept_proposal = accept();

                if (accept_proposal)
                {
                    // store current in the history, replace current by proposal
                    move();
                }
                else
                {
                    // restore previous state of Parameters
                    revert();
                }

                // save points, update statistics etc
                update();
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

        void set_mode(hdf5::File & file, const std::string & data_base_name,
                          const std::vector<double> & point, const double & density)
        {
            // set Stats object
            stats.parameters_at_mode = point;
            stats.mode = density;

            // dump to disk
            auto data_set_mode = file.create_or_open_data_set(data_base_name + "/stats/mode", sample_type);
            std::vector<double> record(parameter_descriptions.size() + 1);
            std::copy(stats.parameters_at_mode.cbegin(), stats.parameters_at_mode.cend(), record.begin());
            record.back() = stats.mode;
            data_set_mode << record;

            return;
        }

        void set_point(const std::vector<double> & point)
        {
            // input validation
            if (parameter_descriptions.size() != point.size())
                throw InternalError("markov_chain::set_point: Dimension of the parameter space of the analysis"
                                    "doesn't match the dimension of the point given.");

            if (parameter_descriptions.empty() || point.empty())
                throw InternalError("markov_chain::set_point: Cannot operate on zero dimensional parameter space");

            {
                auto i = point.cbegin();
                for (auto p = parameter_descriptions.begin(), p_end = parameter_descriptions.end(); p != p_end; ++p, ++i)
                {
                    if (*i < p->min || *i > p->max)
                        throw InternalError("markov_chain::set_point: Parameter '" + p->parameter->name() + "' = " + stringify(*i) + "out of range");
                }
            }

            // assign the values
            {
                for (unsigned i = 0 ; i != parameter_descriptions.size() ; ++i)
                {
                    current.point[i] = point[i];
                    proposal.point[i] = point[i];
                    parameter_descriptions[i].parameter->set(point[i]);
                }
            }

            // copy
            {
                current.log_density = density->evaluate();
                proposal = current;
            }

            // setup statistics
            if (current.log_density > stats.mode)
            {
                stats.mode = current.log_density;
                stats.parameters_at_mode = current.point;
            }

            Log::instance()->message("markov_chain.set_point", ll_debug)
                << current;
        }

        // save points, update statistics
        void update()
        {
            // store points
            if (history.keep)
            {
                history.states.push_back(current);
            }

            if (accept_proposal)
            {
                ++stats.iterations_accepted;
            }
            else
            {
                ++stats.iterations_rejected;
            }

            // count iterations for this parameter since (pre|main) run started. start index at 0, so need +1
            double total_iterations_since_reset = stats.iterations_total + (current_iteration + 1.0);

            if (current.log_density > stats.mode)
            {
                stats.mode = current.log_density;
                stats.parameters_at_mode = current.point;
            }

            // todo remove: let clients figure out means and variances
            for (unsigned i = 0 ; i < parameter_descriptions.size() ; ++i)
            {
                // update mean values, keep copy for variance calculation below

                double former_mean_of_parameter = stats.mean_of_parameters[i];
                stats.mean_of_parameters[i] += (current.point[i] - former_mean_of_parameter) / total_iterations_since_reset;

                if (total_iterations_since_reset < 2)
                {
                    welford_data_parameters[i] = 0;
                }
                else
                {
                    // update variance using Welford's method
                    // see http://www.johndcook.com/standard_deviation.html,
                    // Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd edition
                    welford_data_parameters[i] += (current.point[i] - former_mean_of_parameter) *
                        (current.point[i] - stats.mean_of_parameters[i]);

                    stats.variance_of_parameters[i] = welford_data_parameters[i] / (total_iterations_since_reset - 1);
                }
            }

            // update density
            double former_density = stats.mean_of_log_density;
            stats.mean_of_log_density += (current.log_density - former_density) / total_iterations_since_reset;
            if (total_iterations_since_reset < 2)
            {
                welford_data_density = 0;
            }
            else
            {
                // update variance using Welford's method
                // see http://www.johndcook.com/standard_deviation.html,
                // Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd edition
                welford_data_density += (current.log_density - former_density) * (current.log_density - stats.mean_of_log_density);
                stats.variance_of_log_density = welford_data_density / (total_iterations_since_reset - 1);
            }
        }

        // return random number in [0,1[
        inline double uniform_random_number() { return gsl_rng_uniform(rng); }
    };

    MarkovChain::MarkovChain(const DensityPtr & density, unsigned long seed, const std::shared_ptr<MarkovChain::ProposalFunction> & proposal_function) :
        PrivateImplementationPattern<MarkovChain>(new Implementation<MarkovChain>(density, seed, proposal_function))
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
    MarkovChain::dump_history(hdf5::File & file, const std::string & data_set_base_name, const unsigned & last_iterations) const
    {
        _imp->dump_history(file, data_set_base_name, last_iterations);
    }

    void
    MarkovChain::dump_proposal(hdf5::File & file, const std::string & data_set) const
    {
        _imp->dump_proposal(file, data_set);
    }

    void
    MarkovChain::set_point(const std::vector<double> & point)
    {
        _imp->set_point(point);
    }

    const MarkovChain::State &
    MarkovChain::current_state() const
    {
        return _imp->current;
    }

    const MarkovChain::State &
    MarkovChain::proposed_state() const
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

    ProposalFunctionPtr
    MarkovChain::proposal_function() const
    {
        return _imp->proposal_function;
    }

    void
    MarkovChain::proposal_function(const ProposalFunctionPtr & prop)
    {
        _imp->proposal_function = prop;
    }

    void
    MarkovChain::read_data(hdf5::File & file, const std::string & data_base_name,
                       MarkovChain::History & history, ProposalFunctionPtr & proposal,
                       std::string & proposal_type, MarkovChain::Stats & stats)
    {
        // extract meta information only once
        auto meta_record = proposal_functions::meta_record();
        auto meta_data_set = file.open_data_set(data_base_name + "/proposal/meta", proposal_functions::meta_type());
        meta_data_set >> meta_record;
        proposal_type = std::get<0>(meta_record);
        const unsigned & dimension = std::get<1>(meta_record);

        Implementation<MarkovChain>::read_history(file, data_base_name, dimension, history);
        Implementation<MarkovChain>::read_proposal(file, data_base_name + "/proposal", proposal_type, dimension, proposal);
        Implementation<MarkovChain>::read_stats(file, data_base_name, dimension, stats);
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

    void
    MarkovChain::set_mode(hdf5::File & file, const std::string & data_base_name,
                          const std::vector<double> & point, const double & density)
    {
        _imp->set_mode(file, data_base_name, point, density);
    }

    const MarkovChain::Stats &
    MarkovChain::statistics() const
    {
        return _imp->stats;
    }

    const MarkovChain::History &
    MarkovChain::history() const
    {
        return _imp->history;
    }

    MarkovChain::ProposalFunction::~ProposalFunction()
    {
    }

    const MarkovChain::State &
    MarkovChain::History::local_mode(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end) const
    {
        return *std::max_element(begin, end, [](const MarkovChain::State & a, const MarkovChain::State & b) { return a.log_density < b.log_density; });
    }

    void
    MarkovChain::History::mean_and_variance(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                            std::vector<double> & mean, std::vector<double> & variance) const
    {
        // exclude trivial case
        if (begin == end)
            throw InternalError("MarkovChain::History::mean_and_variance: Cannot compute statistics for empty sequence");

        std::vector<double> temp_means(begin->point.size(), 0.0);

        // initialization to zero
        std::vector<double> temp_squared_sum(begin->point.size(), 0.0);

        // input can be fixed to the right size, and first step calculated in one go
        mean = begin->point;
        variance.assign(begin->point.size(), 0.0);

        // we start at second sample
        unsigned number_of_states = 2;

        // loop over states
        for (auto s = begin + 1 ; s != end ; ++s, ++number_of_states)
        {
            // parameter index
            unsigned i = 0;
            // loop over parameters
            for (auto p = s->point.begin(), p_end = s->point.end() ; p != p_end ; ++p, ++i)
            {
                // calculate the running mean
                temp_means[i] = mean[i];
                mean[i] += (*p - temp_means[i]) / number_of_states;

                // running variance
                temp_squared_sum[i] += (*p - temp_means[i]) * (*p - mean[i]);
                variance[i] = temp_squared_sum[i] / (number_of_states - 1);
            }
        }
    }

    void
    MarkovChain::History::mean_and_covariance(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                            std::vector<double> & mean, std::vector<double> & covariance) const
    {
        const unsigned & dim = begin->point.size();
        std::vector<double> variance(dim, 0.0);

        this->mean_and_variance(begin, end, mean, variance);

        covariance.assign(dim * dim, 0.0);

        // copy diagonal elements
        for (unsigned i = 0 ; i < dim; ++i)
        {
            covariance[i * dim + i] = variance[i];
        }

        const unsigned number_of_history_states = std::distance(begin, end);

        // covariance calculation for off-diagonal elements
        for (auto s = begin ; s != end ; ++s)
        {
            for (unsigned i = 0 ; i < dim; ++i)
            {
                // off-diagonal elements
                for (unsigned j = i + 1 ; j < dim ; ++j)
                {
                    // rescale for the unbiased estimate of sample covariance
                    const double summand = (s->point[i] - mean[i]) * (s->point[j] -  mean[j]);
                    covariance[i + dim * j] += summand;
                    covariance[j + dim * i] += summand;
                }
            }
        }

        // rescale for the unbiased estimate of sample covariance
        for (unsigned i = 0 ; i < dim; ++i)
        {
            for (unsigned j = i + 1 ; j < dim ; ++j)
            {
                covariance[i + dim * j] /= (number_of_history_states - 1);
                covariance[j + dim * i] /= (number_of_history_states - 1);
            }
        }
    }

    std::ostream &
    operator<< (std::ostream & lhs, const MarkovChain::State & rhs)
    {
        lhs << "point = ( ";

        for (auto p = rhs.point.cbegin(), p_end = rhs.point.cend() ; p != p_end ; ++p)
        {
            lhs << *p << " ";
        }

        lhs << "), log(density) = " << rhs.log_density;

        return lhs;
    }
}
