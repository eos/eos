/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <test/test.hh>
#include <src/utils/likelihood.hh>
#include <src/utils/power_of.hh>

#include <cmath>

#include <limits>

using namespace test;
using namespace eos;

struct TestObservable :
    public Observable
{
    Parameters p;

    Kinematics k;

    Options o;

    std::string n, mass_name;

    TestObservable(const Parameters & p, const Kinematics & k, const std::string & mass_name) :
        p(p),
        k(k),
        o(),
        n("test-observable[" + mass_name + "]"),
        mass_name(mass_name)
    {
    }

    double evaluate() const
    {
        double m = p[mass_name];
        double s = k["s"];

        return std::log(power_of<2>(m) / s);
    }

    ObservablePtr clone() const
    {
        static const bool should_never_be_used = true;
        TEST_CHECK(should_never_be_used == false);
    }

    ObservablePtr clone(const Parameters & parameters) const
    {
        return ObservablePtr(new TestObservable(parameters, k.clone(), mass_name));
    }

    Parameters parameters() { return p; }
    Kinematics kinematics() { return k; }
    Options options() { return o; }
    const std::string & name() const { return n; }
};

class LikelihoodTest :
    public TestCase
{
    public:
        LikelihoodTest() :
            TestCase("likelihood_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["mass::b(MSbar)"] = 4.2;
            p["mass::c"] = 1.27;
            p["mass::tau"] = 1.779;

            Kinematics k;
            k.declare("s", 15.0);

            static double eps = std::numeric_limits<double>::epsilon();
            // Zero test
            {
                Likelihood llh(p);
                llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +0.15, +0.162118849476435, +0.17);

                TEST_CHECK_NEARLY_EQUAL(llh(), 1.000000000000000, eps);
            }

            // Single test
            {
                Likelihood llh(p);
                llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +0.15, +0.16, +0.17);

                TEST_CHECK_NEARLY_EQUAL(llh(), 0.994403813437738, eps);
            }

            // Multiple test
            {
                Likelihood llh(p);
                llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +0.15, +0.16, +0.17);
                llh.add(ObservablePtr(new TestObservable(p, k, "mass::c")),        -1.33, -1.82, -1.90);
                llh.add(ObservablePtr(new TestObservable(p, k, "mass::tau")),      -1.85, -2.00, -2.18);

                TEST_CHECK_NEARLY_EQUAL(llh(), 0.310470594301509861, eps);
            }

            // Clone test
            {
                Likelihood llh1(p);
                llh1.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +0.15, +0.16, +0.17);

                TEST_CHECK_NEARLY_EQUAL(llh1(), 0.994403813437738, eps);

                Likelihood llh2 = llh1.clone();
                TEST_CHECK_NEARLY_EQUAL(llh1(), 0.994403813437738, eps);
                TEST_CHECK_NEARLY_EQUAL(llh2(), 0.994403813437738, eps);

                p["mass::b(MSbar)"] = 4.30;
                TEST_CHECK_NEARLY_EQUAL(llh1(), 0.048639401002063827, eps);
                TEST_CHECK_NEARLY_EQUAL(llh2(), 0.994403813437738, eps);

                llh2.parameters()["mass::b(MSbar)"] = 4.30;
                TEST_CHECK_NEARLY_EQUAL(llh1(), 0.048639401002063827, eps);
                TEST_CHECK_NEARLY_EQUAL(llh2(), 0.048639401002063827, eps);
            }
        }
} likelihood_test;
