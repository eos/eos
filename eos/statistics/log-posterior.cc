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

#include <eos/statistics/log-posterior.hh>
#include <eos/utils/density-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnSimplex.h>

#include <gsl/gsl_cdf.h>

using namespace ROOT::Minuit2;

namespace eos
{
    struct RangeError :
        public Exception
    {
        RangeError(const std::string & message) throw () :
            Exception("Range Error: " + message)
        {
        }
    };

   struct MinuitAdapter :
       public ROOT::Minuit2::FCNBase
   {
       LogPosterior log_posterior;

       MnUserParameters user_parameters;

       // stores results of minimization
       std::shared_ptr<FunctionMinimum> data_at_minimum;

       MinuitAdapter(LogPosterior & log_posterior) :
           log_posterior(log_posterior)
       {
           // define parameters and limits
           for (auto i = log_posterior.parameter_descriptions().begin(), i_end = log_posterior.parameter_descriptions().end() ; i != i_end ; ++i)
           {
               // name, value, error(?), min, max
               user_parameters.Add(i->parameter->name(), (i->min + i->max) / 2.0,  i->max - i->min, i->min, i->max);
           }
       }

       virtual ~MinuitAdapter()
       {
       }

       virtual double Up() const
       {
           return 0.5;
       }

       virtual double operator()(const std::vector<double> & parameter_values) const
       {
           // copy parameter values
           auto p = parameter_values.begin();
           for (auto i = log_posterior.parameter_descriptions().begin(), i_end = log_posterior.parameter_descriptions().end() ; i != i_end ; ++i, ++p)
           {
               i->parameter->set(*p);
           }
           return -log_posterior.log_posterior();
       }
   };

   LogPosterior::LogPosterior(const LogLikelihood & log_likelihood) :
        _log_likelihood(log_likelihood),
        _parameters(log_likelihood.parameters()),
        _minuit(nullptr)
   {
   }

   LogPosterior::~LogPosterior()
   {
       delete _minuit;
   }

   bool
   LogPosterior::add(const LogPriorPtr & prior, bool nuisance)
   {
       // clone has correct Parameters object selected
       LogPriorPtr prior_clone = prior->clone(_parameters);

       // check if param exists already
       // read out parameters from prior
       for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d)
       {
           auto result = _parameter_names.insert(d->parameter->name());
           if (! result.second)
               return false;

           d->nuisance = nuisance;
           _parameter_descriptions.push_back(*d);
       }

       // then add to prior container
       _priors.push_back(prior_clone);

       return true;
   }

   DensityPtr
   LogPosterior::clone() const
   {
       return DensityPtr(private_clone());
   }

   LogPosteriorPtr
   LogPosterior::old_clone() const
   {
       return LogPosteriorPtr(private_clone());
   }

   LogPosterior *
   LogPosterior::private_clone() const
   {
       // clone log_likelihood
       LogLikelihood llh = _log_likelihood.clone();
       LogPosterior * result = new LogPosterior(llh);

       // add parameters via prior clones
       for (auto i = _priors.cbegin(), i_end = _priors.cend(); i != i_end; ++i)
       {
           result->add((*i)->clone(result->parameters()));
       }

       // copy proper range for subspace sampling
       auto j = result->_parameter_descriptions.begin();
       for (auto i = _parameter_descriptions.begin(), i_end = _parameter_descriptions.end() ; i != i_end ; ++i, ++j)
       {
           j->min = i->min;
           j->max = i->max;
           j->nuisance = i->nuisance;
       }

       return result;
   }

   double
   LogPosterior::evaluate() const
   {
       return log_posterior();
   }

   Density::Iterator
   LogPosterior::begin() const
   {
       return Density::Iterator(_parameter_descriptions.cbegin());
   }

   Density::Iterator
   LogPosterior::end() const
   {
       return Density::Iterator(_parameter_descriptions.cend());
   }

   void
   LogPosterior::dump_descriptions(hdf5::File & file, const std::string & data_set_base) const
   {
       // store parameter info, including the prior
       {
           auto data_set = file.create_data_set(data_set_base + "/parameters", LogPosterior::Output::description_type());

           auto record = LogPosterior::Output::description_record();
           std::string prior;

           for (auto d = _parameter_descriptions.cbegin(), d_end = _parameter_descriptions.cend() ; d != d_end ; ++d)
           {
               std::get<0>(record) = d->parameter->name().c_str();
               std::get<1>(record) = d->min;
               std::get<2>(record) = d->max;
               std::get<3>(record) = int(d->nuisance);

               // without local string variable, can get random data in char *
               prior = log_prior(d->parameter->name())->as_string();
               std::get<4>(record) = prior.c_str();

               data_set << record;
           }
           // store the SHA hash of the current git version
           auto attr = data_set.create_attribute("version", hdf5::Scalar<const char *>("version"));
           attr = EOS_GITHEAD;
       }

       // store constraints
       {
           hdf5::Composite<hdf5::Scalar<const char *>> constraint_type
           {
               "constraints",
               hdf5::Scalar<const char *>("name"),
           };
           auto constraint_data_set = file.create_data_set(data_set_base + "/constraints", constraint_type);
           auto constraint_record = std::make_tuple("name");

           // loop over constraints
           for (const auto & c : _log_likelihood)
           {
               std::get<0>(constraint_record) = c.name().str().c_str();
               constraint_data_set << constraint_record;
           }
       }

       // store observable names
       {
           hdf5::Composite<hdf5::Scalar<const char *>> observables_type
           {
               "observables",
               hdf5::Scalar<const char *>("name"),
           };

           auto observables_data_set = file.create_data_set(data_set_base + "/observables", observables_type);
           auto observables_record = std::make_tuple("name");

           const ObservableCache & cache = _log_likelihood.observable_cache();
           const unsigned & n_observables = cache.size();
           for (unsigned i = 0 ; i < n_observables ; ++i)
           {
               std::get<0>(observables_record) = cache.observable(i)->name().str().c_str();
               observables_data_set << observables_record;
           }
       }
   }

   std::vector<ParameterDescription>
   LogPosterior::read_descriptions(const hdf5::File & file, std::string data_set_base)
   {
       hdf5::File & f = const_cast<hdf5::File &>(file);
       auto data_set = f.open_data_set(data_set_base + "/parameters", LogPosterior::Output::description_type());

       auto record = LogPosterior::Output::description_record();

       std::vector<ParameterDescription> descriptions;
       Parameters p = Parameters::Defaults();

       for (unsigned i = 0 ; i < data_set.records() ; ++i)
       {
           data_set >> record;

           descriptions.push_back(ParameterDescription { p[std::get<0>(record)].clone(), std::get<1>(record), std::get<2>(record), bool(std::get<3>(record)) });
       }
       return descriptions;
   }

   void
   LogPosterior::read_descriptions(hdf5::File & file, const std::string & data_set_base_name,
           std::vector<ParameterDescription> & descr, std::vector<std::string> & priors,
           std::vector<std::string> & constraints, std::string & hash)
   {
       hdf5::Composite<hdf5::Scalar<const char *>, hdf5::Scalar<double>,
       hdf5::Scalar<double>, hdf5::Scalar<int>, hdf5::Scalar<const char *>> parameter_descriptions_type
       {
           "parameter description",
           hdf5::Scalar<const char *>("name"),
           hdf5::Scalar<double>("min"),
           hdf5::Scalar<double>("max"),
           hdf5::Scalar<int>("nuisance"),
           hdf5::Scalar<const char *>("prior"),
       };
       auto data_set = file.open_data_set(data_set_base_name + "/parameters", parameter_descriptions_type);

       Parameters p = Parameters::Defaults();

       for (unsigned i = 0 ; i < data_set.records() ; ++i)
       {
           auto record = std::make_tuple("parameter_name", 1.0, 2.0, 3, "prior");
           data_set >> record;
           descr.push_back(ParameterDescription{ p[std::get<0>(record)].clone(), std::get<1>(record),
               std::get<2>(record), bool(std::get<3>(record)) });
           priors.push_back(std::get<4>(record));
       }

       // store the SHA hash of the current git version
       auto attr = data_set.open_attribute("version", hdf5::Scalar<const char *>("version"));
       hash = attr.value();

       // store constraints
       {
           hdf5::Composite<hdf5::Scalar<const char *>> constraint_type
           {
               "constraints",
               hdf5::Scalar<const char *>("name"),
           };
           auto constraint_data_set = file.open_data_set(data_set_base_name + "/constraints", constraint_type);
           auto constraint_record = std::make_tuple("name");

           for (unsigned i = 0 ; i < constraint_data_set.records() ; ++i)
           {
               constraint_data_set >> constraint_record;
               constraints.push_back(std::get<0>(constraint_record));
           }
       }
   }

   Parameters
   LogPosterior::parameters() const
   {
       return _parameters;
   }

   std::pair<double, double>
   LogPosterior::goodness_of_fit(const std::vector<double> & parameter_values, const unsigned & simulated_datasets, std::string output_file)
   {
       // count scan parameters
       unsigned scan_parameters = 0;

       for (auto d = _parameter_descriptions.cbegin() ; d != _parameter_descriptions.cend() ; ++d)
       {
           if (! d->nuisance)
               ++scan_parameters;
       }

       if (_parameter_descriptions.size() != parameter_values.size())
           throw InternalError("LogPosterior::goodness_of_fit: starting point doesn't have the correct dimension: "
                   + stringify(parameter_values.size()) + " vs " + stringify(_parameter_descriptions.size()) );

       std::shared_ptr<hdf5::File> f;
       if ( ! output_file.empty())
       {
           f = std::shared_ptr<hdf5::File>(new hdf5::File(hdf5::File::Create(output_file)));
           dump_descriptions(*f, "/descriptions");

           auto data_set = f->create_data_set("/data/parameters", hdf5::Array<1, double>("goodness-of-fit-point", { parameter_values.size() }));
           data_set << parameter_values;
       }

       // set the parameter values
       for (unsigned i = 0 ; i < parameter_values.size(); ++i)
       {
           ParameterDescription & desc = _parameter_descriptions[i];
           if (parameter_values[i] < desc.min || parameter_values[i] > desc.max)
           {
               throw InternalError("LogPosterior::goodness_of_fit: parameter " + desc.parameter->name() +
                       " out of bounds [" + stringify(desc.min) + ", " + stringify(desc.max) + "]: "
                       + stringify(parameter_values[i]));
           }
           desc.parameter->set(parameter_values[i]);
       }

       // update observables for new parameter values
       const double ll = _log_likelihood();

       Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                        << "Calculating p-values at parameters " << stringify(parameter_values.cbegin(), parameter_values.cend())
                        << " with log(post) = " << ll + log_prior();

       // simulate data sets
       auto sim_result = _log_likelihood.bootstrap_p_value(simulated_datasets);

       // p-value from the analytical, yet approximate \chi^2-distribution
       // with (n_obs - n_par) degrees-of-freedom
       const double dof = 1.0 * _log_likelihood.number_of_observations() - _parameter_descriptions.size();
       const double chi_squared = gsl_cdf_chisq_Qinv(sim_result.first, _log_likelihood.number_of_observations());

       double p_analytical = 0;
       if (dof > 0)
       {
           p_analytical = gsl_cdf_chisq_Q(chi_squared, dof);

           Log::instance()->message("log_posterior.goodness_of_fit", ll_debug)
                                << "dof = " << dof << ", _parameter_descriptions.size = " << _parameter_descriptions.size()
                                << ", #observations = " << _log_likelihood.number_of_observations();

           Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                               << "p-value from simulating pseudo experiments after applying DoF correction and using the \\chi^2-distribution"
                               << " (valid assumption?) has a value of " << p_analytical;
       }
       else
       {
           Log::instance()->message("log_posterior.goodness_of_fit", ll_warning)
                            << "Cannot compute p-value for non-positive dof (" << dof << "). Need more constraints / less parameters";
       }

       // p-value from the analytical, yet approximate \chi^2-distribution
       // with (n_obs - n_par) degrees-of-freedom
       double dof_scan = 1.0 * _log_likelihood.number_of_observations() - scan_parameters;
       if (dof_scan > 0)
       {
           double p_analytical_scan = gsl_cdf_chisq_Q(chi_squared, dof_scan);

           Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                           << "p-value from simulating pseudo experiments after applying DoF correction (scan parameters only)"
                           << " and using the \\chi^2-distribution"
                           << " (valid assumption?) has a value of " << p_analytical_scan;
       }
       else
       {
           Log::instance()->message("log_posterior.goodness_of_fit", ll_warning)
                            << "Cannot compute p-value for negative dof_scan (" << dof_scan << "). Need more constraints / less parameters";

       }
       /* calculate the significances */
       double total_significance_squared = 0;
       std::vector<double> significances;

       Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                        << "Significances for each constraint:";

       for (auto c = _log_likelihood.begin(), c_end = _log_likelihood.end() ; c != c_end ; ++c)
       {
           for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
           {
               double significance = (**b).significance();
               Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                                << c->name() << ": " << significance << " sigma";
               total_significance_squared += power_of<2>(significance);
               significances.push_back(significance);
           }
       }

       /* output the test statistics */
       Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                        << "TestStatistic for each constraint:";
       for (auto c = _log_likelihood.begin(), c_end = _log_likelihood.end() ; c != c_end ; ++c)
       {
           for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
           {
               std::stringstream ss;
               ss << *(*b)->primary_test_statistic();

               Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                                << c->name() << ": " << ss.str();
           }
       }

       // store significances and chi^2
       if (f)
       {
           auto dtype = hdf5::Array<1, double>("goodness-of-fit-point", { significances.size() });
           auto data_set = f->create_data_set("/data/significances", dtype);
           data_set << significances;
           auto attr_sign = data_set.create_attribute("chi2_significance", hdf5::Scalar<double>("chi2_significance"));
           attr_sign = total_significance_squared;

           auto attr_sim = data_set.create_attribute("chi2_simulation", hdf5::Scalar<double>("chi2_simulation"));
           attr_sign = chi_squared;
       }

       Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                        << "Listing the individual observables' predicted values:";

       const ObservableCache & cache = _log_likelihood.observable_cache();
       for (unsigned i = 0 ; i < cache.size() ; ++i)
       {
           Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                            << cache.observable(i)->name() << " = " << cache[i];
       }

       if (dof > 0)
       {
           const double p_significance =  gsl_cdf_chisq_Q(total_significance_squared, dof);
           Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                            << "p-value from calculating significances, treating them as coming from a Gaussian, is "
                            << p_significance << ". The pseudo chi_squared/dof is " << total_significance_squared
                            << "/" << dof << " = " << total_significance_squared / dof;
       }

       if (dof_scan > 0)
       {
           const double p_significance_scan =  gsl_cdf_chisq_Q(total_significance_squared, dof_scan);
           Log::instance()->message("log_posterior.goodness_of_fit", ll_informational)
                                << "p-value from calculating significances, treating them as coming from a Gaussian, is "
                                << p_significance_scan << ". The pseudo chi_squared/dof (dof from scan parameters only) is "
                                << total_significance_squared
                                << "/" << dof_scan << " = " << total_significance_squared / dof_scan;
       }

       return std::make_pair(sim_result.first, p_analytical);
   }

   unsigned
   LogPosterior::index(const std::string & name) const
   {
       unsigned result = 0;

       for (auto d = _parameter_descriptions.cbegin(), d_end = _parameter_descriptions.cend() ; d != d_end ; ++d, ++result)
       {
           if (name == d->parameter->name())
               return result;
       }

       throw InternalError("Implementation<Analysis>::definition: no such parameter '" + name + "'");
   }

   LogLikelihood
   LogPosterior::log_likelihood() const
   {
       return _log_likelihood;
   }

   double
   LogPosterior::log_posterior() const
   {
       return log_prior() + _log_likelihood();
   }

   double
   LogPosterior::log_prior() const
   {
       if (_priors.empty())
           throw InternalError("LogPosterior::log_prior(): prior is undefined");

       double result = 0.0;

       // all prior components are assumed independent,
       // thus the logs can be simply added up
       for (auto p = _priors.cbegin(), p_end = _priors.cend() ; p != p_end; ++p)
       {
           result += (**p)();
       }

       return result;
   }

   LogPriorPtr
   LogPosterior::log_prior(const std::string & name) const
   {
       LogPriorPtr prior;

       // loop over all descriptions of the prior pointers
       for (auto p = _priors.begin(), p_end = _priors.end() ; p != p_end ; ++p)
       {
           for (auto i = (*p)->begin(), i_end = (*p)->end() ; i != i_end ; ++i)
           {
               if (i->parameter->name() == name)
                   prior = *p;
           }
       }

       return prior;
   }

   double
   LogPosterior::negative_log_posterior(const gsl_vector * pars, void * data)
   {
       // set all components of parameters
       LogPosterior * log_posterior = static_cast<LogPosterior *>(data);
       for (unsigned i = 0 ; i < log_posterior->_parameter_descriptions.size() ; ++i)
       {
           log_posterior->_parameter_descriptions[i].parameter->set(gsl_vector_get(pars, i));
       }

       // calculate negative posterior
       return -(log_posterior->log_prior() + log_posterior->log_likelihood()());
   }

   bool
   LogPosterior::nuisance(const std::string& par_name) const
   {
       unsigned index = this->index(par_name);

       if (index >= _parameter_descriptions.size())
       {
           return false;
       }
       else
       {
           return _parameter_descriptions[index].nuisance;
       }
   }

   MutablePtr
   LogPosterior::operator[] (const unsigned & index) const
   {
       return _parameter_descriptions[index].parameter;
   }

   std::pair<std::vector<double>, double>
   LogPosterior::optimize(const std::vector<double> & initial_guess, const LogPosterior::OptimizationOptions & options)
   {
       // input validation
       if (_parameter_descriptions.size() != initial_guess.size())
           throw InternalError("LogPosterior::optimize: starting point doesn't have the correct dimension "
                   + stringify(_parameter_descriptions.size()) );

       // setup the function object as in GSL manual 36.4
       gsl_multimin_function posterior;
       posterior.n = _parameter_descriptions.size();
       posterior.f = &LogPosterior::negative_log_posterior;
       posterior.params = (void *) this;

       // set starting guess
       gsl_vector *x = gsl_vector_alloc(posterior.n);
       for (unsigned i = 0 ; i < posterior.n ; ++i)
           gsl_vector_set(x, i, initial_guess[i]);

       // save minimum for later comparison
       double initial_minimum = posterior.f(x, (void *) this);

       // set initial step sizes relative to allowed parameter range
       gsl_vector *ss = gsl_vector_alloc(posterior.n);
       for (unsigned i = 0 ; i < posterior.n ; ++i)
           gsl_vector_set(ss, i, (_parameter_descriptions[i].max - _parameter_descriptions[i].min) * options.initial_step_size);

       // setup the minimizer
       const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
       gsl_multimin_fminimizer *minim = gsl_multimin_fminimizer_alloc(T, posterior.n);
       gsl_multimin_fminimizer_set(minim, &posterior, x, ss);

       unsigned iter = 0;
       int status = 0;
       double simplex_size = 0;

       // run the minimizer
       do
       {
           iter++;
           status = gsl_multimin_fminimizer_iterate(minim);

           if (status)
               break;

           simplex_size = gsl_multimin_fminimizer_size(minim);
           status = gsl_multimin_test_size(simplex_size, options.tolerance);

           Log::instance()->message("log_posterior.optimize", ll_debug)
                                       << "f() = " << minim->fval << "\tsize = " << simplex_size;

           if (status == GSL_SUCCESS)
           {
               Log::instance()->message("log_posterior.optimize", ll_informational)
                               << "Simplex algorithm converged after " << stringify(iter) << " iterations";
           }
       } while (status == GSL_CONTINUE && iter < options.maximum_iterations);

       // build output vector
       std::vector<double> parameters_at_mode(initial_guess);
       double mode = minim->fval;
       for (unsigned i = 0 ; i < posterior.n ; ++i)
           parameters_at_mode[i] = gsl_vector_get(minim->x, i);

       // free resources
       gsl_vector_free(x);
       gsl_vector_free(ss);
       gsl_multimin_fminimizer_free(minim);

       //check if algorithm actually found a better minimum
       if (mode >= initial_minimum)
       {
           Log::instance()->message("log_posterior.optimize", ll_warning)
                          << "Simplex algorithm did not improve on initial guess";
           return std::make_pair(initial_guess, -initial_minimum) ;
       }

       std::string results("Results: maximum of posterior = ");
       results += stringify(-mode) + " at ( ";
       for (unsigned i = 0 ; i < posterior.n ; ++i)
           results += stringify(parameters_at_mode[i]) + " ";
       results += ")";
       Log::instance()->message("log_posterior.optimize", ll_informational)
                       << results;

       //minus sign to convert to posterior
       return std::make_pair(parameters_at_mode, -mode) ;

   }

   const ROOT::Minuit2::FunctionMinimum &
   LogPosterior::optimize_minuit(const std::vector<double> & initial_guess, const LogPosterior::OptimizationOptions & options)
   {
       if (! _minuit)
       {
           _minuit = new MinuitAdapter(*this);
       }
       // copy values
       MnUserParameters & mn_par(_minuit->user_parameters);
       unsigned j = 0;
       for (auto i = initial_guess.begin(), i_end = initial_guess.end() ; i != i_end ; ++i, ++j)
       {
           mn_par.SetValue(j, *i);
       }

       // fix nuisance parameters with a flat prior to avoid flat directions causing Minuit to fail
       if (options.fix_flat_nuisance)
       {
           // loop over parameters and find those nuisance parameters with flat prior
           unsigned i = 0;
           for (auto d = _parameter_descriptions.cbegin(), i_end = _parameter_descriptions.cend() ; d != i_end ; ++d, ++i)
           {
               if ( ! d->nuisance)
                   continue;

               auto p = log_prior(d->parameter->name());

               // dirty hack: check for flat prior from the variance: works only if this dimension is not partitioned
               if (std::fabs(p->variance() - power_of<2>(d->max - d->min) / 12.0) < 1e-15)
               {
                   // to simplify mode identification, the fixed values agree between different chains
                   mn_par.SetValue(i, d->min);
                   mn_par.Fix(i);
               }
           }
       }

       // create MIGRAD minimizer
       std::shared_ptr<ROOT::Minuit2::MnApplication> minimizer;

       if (options.algorithm == "migrad")
       {
           minimizer.reset(new MnMigrad(*_minuit, mn_par, options.strategy_level));
       }
       // uses minuit, reverts to simplex if it failed, then call minuit again
       if (options.algorithm == "minimize")
       {
           minimizer.reset(new MnMinimize(*_minuit, mn_par, options.strategy_level));
       }
       if (options.algorithm == "scan")
       {
           minimizer.reset(new MnScan(*_minuit, mn_par, options.strategy_level));
       }
       if (options.algorithm == "simplex")
       {
           minimizer.reset(new MnSimplex(*_minuit, mn_par, options.strategy_level));
       }

       if (! minimizer)
           throw InternalError("LogPosterior::optimize_minut: invalid algorithm option: " + options.algorithm);

       // minimize and save results
       _minuit->data_at_minimum.reset(new FunctionMinimum((*minimizer)(options.maximum_iterations, options.tolerance)));

       // release parameters again
       if (options.fix_flat_nuisance)
       {
           unsigned i = 0;
           for (auto d = _parameter_descriptions.cbegin(), i_end = _parameter_descriptions.cend() ; d != i_end ; ++d, ++i)
           {
               mn_par.Release(i);
           }
       }

       return *_minuit->data_at_minimum;
   }

    const std::vector<ParameterDescription> &
    LogPosterior::parameter_descriptions() const
    {
        return _parameter_descriptions;
    }

    LogPosterior::OptimizationOptions::OptimizationOptions():
        algorithm("minimize"),
        fix_flat_nuisance(false),
        initial_step_size(0, 1, 0.1),
        maximum_iterations(8000),
        mcmc_pre_run(true),
        tolerance(0, 1, 1e-1),
        splitting_tolerance(0, 1, 1e-2),
        strategy_level(0,2, 1)
    {
    }

    LogPosterior::OptimizationOptions
    LogPosterior::OptimizationOptions::Defaults()
    {
        return OptimizationOptions();
    }

    LogPosterior::Output::DescriptionType
    LogPosterior::Output::description_type()
    {
        return
        DescriptionType
        {
         "parameter description",
         hdf5::Scalar<const char *>("name"),
         hdf5::Scalar<double>("min"),
         hdf5::Scalar<double>("max"),
         hdf5::Scalar<int>("nuisance"),
         hdf5::Scalar<const char *>("prior"),
        };
    }

    std::tuple<const char *, double, double, int, const char *>
    LogPosterior::Output::description_record()
    {
        return std::make_tuple("name", 1.0, 2.0, 3, "prior");
    }

    std::vector<double>
    proposal_covariance(const LogPosterior & log_posterior,
                        double scale_reduction,
                        bool scale_nuisance)
    {
        const auto & npar = log_posterior.parameter_descriptions().size();

        // zero off-diagonal
        std::vector<double> covariance (npar * npar, 0);

        // prior variance on the diagonal
        unsigned par = 0;
        for (auto & def : log_posterior.parameter_descriptions())
        {
            LogPriorPtr prior = log_posterior.log_prior(def.parameter->name());
            covariance[par + npar * par] = prior->variance();

            // rescale variance of scan parameters with a configurable value, in order
            // to avoid drawing too many samples outside the allowed range.
            if (! def.nuisance || scale_nuisance)
            {
                covariance[par + npar * par] /= power_of<2>(scale_reduction);
            }
            ++par;
        }

        return covariance;
    }
}
