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

#include <eos/utils/log_likelihood_C.hh>

#include <cstring>

#define errorhandler(add_constraint_call) \
    do \
    { \
        std::string s(""); \
        try \
        { \
            add_constraint_call; \
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

    LogLikelihood *
    EOS_LogLikelihood_new()
    {
        return new LogLikelihood(Parameters::Defaults());
    }

    void
    EOS_LogLikelihood_delete(LogLikelihood * ll)
    {
        delete ll;
        ll = nullptr;
    }

    char *
    EOS_LogLikelihood_add_constraint_by_name(LogLikelihood * ll,
                                             const char * constraint_name,
                                             const Options * options)
    {
        errorhandler(ll->add(Constraint::make(constraint_name, *options)));
    }

    char *
    EOS_LogLikelihood_add_gaussian_constraint(LogLikelihood * ll, const char * observable,
                                              const double min, const double central,
                                              const double max, const unsigned number_of_observations,
                                              const Kinematics * kinematics, const Options * options)
    {
        errorhandler(
                // create observable
                auto o = Observable::make(observable, ll->parameters(), *kinematics, *options);
                if (not o)
                    throw ObservableNameError(std::string("Unknown observable: ") + observable);
                ll->add(o, min, central, max, number_of_observations);
                );
    }
}
