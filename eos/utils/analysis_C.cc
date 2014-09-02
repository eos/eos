/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Stephan Jahn
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

#include <eos/utils/analysis_C.hh>
#include <eos/utils/log.hh>

#include <cstring>
#include <iostream>
#include <sstream>

#define errorhandler(prior_constructor_call) \
    do \
    { \
        std::string s(""); \
        try \
        { \
            prior_constructor_call; \
            bool added = ana->add(prior, nuisance); \
            if (!added) \
            { \
                s = std::string("EOS: Attempting to add parameter \"") + par_name + "\" twice."; \
            } \
        } \
        catch (eos::Exception & e) \
        { \
            s = "EOS: "; \
            s += e.what(); \
        } \
        catch (...) \
        { \
            s = "Unknown Error"; \
        } \
        /* add "sizeof(char)" to "s.size()" to make sure that there is memory for "NULL" at end of string \
         * Note: A std::string is NOT neccessarily NULL-terminated \
         */ \
        char * c = static_cast<char *>(malloc(s.size() + sizeof(char))); \
        strcpy(c, s.c_str()); \
        return c; \
    } \
    while (false)

extern "C" {
    using namespace eos;

    Analysis *
    EOS_Analysis_new(LogLikelihood * log_likelihood)
    {
        return new Analysis(*log_likelihood);
    }

    void
    EOS_Analysis_delete(Analysis * ana)
    {
        delete ana;
        ana = nullptr;
    }

    char *
    EOS_Analysis_add_Flat(Analysis * ana, const char * par_name,
                          const double range_min, const double range_max,
                          bool nuisance)
    {
        errorhandler(LogPriorPtr prior = LogPrior::Flat(ana->log_likelihood().parameters(), par_name,
                                                        ParameterRange{range_min, range_max}));
    }

    char *
    EOS_Analysis_add_Gauss(Analysis * ana, const char * par_name,
                           const double range_min, const double range_max,
                           const double lower, const double central, const double upper,
                           bool nuisance)
    {
        errorhandler(LogPriorPtr prior = LogPrior::Gauss(ana->log_likelihood().parameters(), par_name,
                                                         ParameterRange{range_min, range_max},
                                                         lower, central, upper) );
    }

    char *
    EOS_Analysis_add_LogGamma(Analysis * ana, const char * par_name,
                              const double range_min, const double range_max,
                              const double lower, const double central, const double upper,
                              bool nuisance)
    {
        errorhandler(LogPriorPtr prior = LogPrior::LogGamma(ana->log_likelihood().parameters(), par_name,
                                                            ParameterRange{range_min, range_max}, lower, central, upper) );
    }

    char *
    EOS_Analysis_info(Analysis * ana)
    {
        std::stringstream ostream;

        ostream << "Constraints:" << std::endl;
        ostream << "------------" << std::endl;
        auto likelih = ana->log_likelihood();
        for (auto & constr : likelih)
            ostream << constr.name() << std::endl;
        ostream << std::endl;

        ostream << std::endl << "Observables:";
        ostream << std::endl << "------------" << std::endl;
        const auto & cache = likelih.observable_cache();
        for (ObservableCache::Id i = 0; i < cache.size(); ++i)
        {
            auto o = cache.observable(i);
            ostream << o->name();
            ostream << ( o->kinematics().as_string() == "" ? "" : ( "[" + o->kinematics().as_string() + "]" )  );
            ostream << " = " << cache[i];
            ostream << std::endl;
            ostream << "  options: " << o->options().as_string();
            ostream << std::endl << std::endl;
        }
        ostream << std::endl;

        ostream << "Parameters:" << std::endl;
        ostream << "-----------" << std::endl;
        for (auto & i : ana->parameter_descriptions())
        {
            Parameter p = i.parameter;
            ostream << ana->log_prior(p.name())->as_string();
            ostream << ", value = " << p() << std::endl << std::endl;
        }
        ostream << std::endl;

        std::string s = ostream.str();

        /* add "sizeof(char)" to "s.size()" to make sure that there is memory for "NULL" at end of string
         * Note: A std::string is NOT neccessarily NULL-terminated
         */
        char * c = (char *) malloc(s.size() + sizeof(char));
        strcpy(c, s.c_str());
        return c;
    }

    char *
    EOS_Analysis_gof(Analysis * ana, const double * par_vals, const unsigned simulated_datasets)
    {
        Log::instance()->set_program_name("eos.py");
        Log::instance()->set_log_level(ll_debug);

        // catch the output
        std::stringstream ostream;
        Log::instance()->set_log_stream(&ostream);

        // copy parameter values into vector
        std::vector<double> par_vals_(par_vals, par_vals + ana->parameter_descriptions().size());

        // ignore return value
        ana->goodness_of_fit(par_vals_, simulated_datasets, "");

        std::string s = ostream.str();

        /* add "sizeof(char)" to "s.size()" to make sure that there is memory for "NULL" at end of string
         * Note: A std::string is NOT neccessarily NULL-terminated
         */
        char * c = (char *) malloc(s.size() + sizeof(char));
        strcpy(c, s.c_str());
        return c;
    }

    double
    EOS_Analysis_log_posterior(Analysis* ana, const double * par_vals)
    {
        unsigned j = 0;
        for (auto & i : ana->parameter_descriptions())
        {
            Parameter p = i.parameter;
            p = par_vals[j];
            ++j;
        }
        // evaluate likelihood and prior; no NaN check --> to be done in python
        return ana->log_posterior();
    }
}
