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

#ifndef EOS_GUARD_SRC_UTILS_ANALYSIS_TEST_HH
#define EOS_GUARD_SRC_UTILS_ANALYSIS_TEST_HH 1

#include <src/utils/analysis.hh>
#include <test/test.hh>

using namespace test;

namespace eos
{
    struct TestObservable :
        public Observable
    {
            Parameters p;

            Kinematics k;

            Options o;

            std::string n, mass_name;

            UsedParameter mass;

            TestObservable(const Parameters & p, const Kinematics & k, const std::string & mass_name) :
                p(p),
                k(k),
                o(),
                n("test-observable[" + mass_name + "]"),
                mass_name(mass_name),
                mass(p[mass_name], *this)
            {
            }

            double evaluate() const
            {
                return mass();
            }

            ObservablePtr clone() const
            {
                return ObservablePtr(new TestObservable(p.clone(), k.clone(), mass_name));
            }

            ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new TestObservable(parameters, k.clone(), mass_name));
            }

            Parameters parameters()
            {
                return p;
            }

            Kinematics kinematics()
            {
                return k;
            }

            Options options()
            {
                return o;
            }

            const std::string & name() const
            {
                return n;
            }

            void set_option(const std::string & key, const std::string & value = "")
            {
                o.set(key, value);
            }
    };

    /*
     * Create analysis with gaussian likelihood and gaussian prior
     * The posterior is also a Gaussian with central value 4.3
     * and standard deviation sqrt(0.005)=0.070710678118654752
     *
     * flat = true turns Gaussian prior into flat prior
     */
    Analysis make_analysis(bool flat)
    {
        Parameters parameters = Parameters::Defaults();

        LogLikelihoodPtr llh(new LogLikelihood(parameters));
        llh->add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                        "mass::b(MSbar)")), 4.1, 4.2, 4.3);

        LogPriorPtr prior = flat ?
            LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{3.7, 4.9} ) :
            LogPrior::Gauss(parameters, "mass::b(MSbar)", ParameterRange{3.7, 4.9}, 4.3, 4.4, 4.5);

        Analysis result(llh);
        result.add(prior);

        return result;
    }
}

#endif
