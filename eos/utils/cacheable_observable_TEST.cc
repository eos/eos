/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2021-2026 Danny van Dyk
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

#include "eos/utils/kinematic.hh"
#include <eos/observable.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

namespace eos
{

    class TestCacheableObservableProvider : public ParameterUser, public PrivateImplementationPattern<TestCacheableObservableProvider>
    {
        public:
            struct IntermediateResult : public CacheableObservable::IntermediateResult
            {
                    // e.g. amplitudes
                    double a;
                    double b;

                    // e.g. kinematical variables
                    double q2;
            };

            struct AlternativeResult : public CacheableObservable::IntermediateResult
            {
                    double c;
            };

            TestCacheableObservableProvider(const Parameters & parameters, const Options & options);
            ~TestCacheableObservableProvider();

            // Observables
            const IntermediateResult * prepare(const double & q2) const;
            const IntermediateResult * prepare_shifted(const double & q2) const;
            const AlternativeResult *  prepare_alternative(const double & q2) const;

            double evaluate1(const IntermediateResult *) const;
            double evaluate2(const IntermediateResult *) const;
            double evaluate3(const AlternativeResult *) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };

    const std::set<ReferenceName> TestCacheableObservableProvider::references{};

    template <> struct Implementation<TestCacheableObservableProvider>
    {
            UsedParameter m_B;

            using IntermediateResult = TestCacheableObservableProvider::IntermediateResult;
            using AlternativeResult  = TestCacheableObservableProvider::AlternativeResult;

            IntermediateResult _intermediate_result;

            AlternativeResult _alternative_result;

            Implementation(const Parameters & p, const Options & /* o */, ParameterUser & u) :
                m_B(p["mass::B_u"], u)
            {
            }

            const IntermediateResult *
            prepare(const double & q2)
            {
                _intermediate_result.a = 2.0;
                _intermediate_result.b = m_B;

                _intermediate_result.q2 = q2;

                return &_intermediate_result;
            }

            // shares its result type and its arity with prepare(), but is a different function
            const IntermediateResult *
            prepare_shifted(const double & q2)
            {
                return prepare(q2 + 1.0);
            }

            const AlternativeResult *
            prepare_alternative(const double & q2)
            {
                _alternative_result.c = m_B * q2;

                return &_alternative_result;
            }

            double
            evaluate1(const IntermediateResult * intermediate_result)
            {
                return intermediate_result->b - intermediate_result->a * intermediate_result->q2;
            }

            double
            evaluate2(const IntermediateResult * intermediate_result)
            {
                return pow(intermediate_result->q2, 2);
            }

            double
            evaluate3(const AlternativeResult * alternative_result)
            {
                return alternative_result->c;
            }
    };

    TestCacheableObservableProvider::TestCacheableObservableProvider(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TestCacheableObservableProvider>(new Implementation<TestCacheableObservableProvider>(parameters, options, *this))
    {
    }

    TestCacheableObservableProvider::~TestCacheableObservableProvider() {}

    const TestCacheableObservableProvider::IntermediateResult *
    TestCacheableObservableProvider::prepare(const double & q2) const
    {
        return _imp->prepare(q2);
    }

    const TestCacheableObservableProvider::IntermediateResult *
    TestCacheableObservableProvider::prepare_shifted(const double & q2) const
    {
        return _imp->prepare_shifted(q2);
    }

    const TestCacheableObservableProvider::AlternativeResult *
    TestCacheableObservableProvider::prepare_alternative(const double & q2) const
    {
        return _imp->prepare_alternative(q2);
    }

    double
    TestCacheableObservableProvider::evaluate1(const TestCacheableObservableProvider::IntermediateResult * ir) const
    {
        return _imp->evaluate1(ir);
    }

    double
    TestCacheableObservableProvider::evaluate2(const TestCacheableObservableProvider::IntermediateResult * ir) const
    {
        return _imp->evaluate2(ir);
    }

    double
    TestCacheableObservableProvider::evaluate3(const TestCacheableObservableProvider::AlternativeResult * ir) const
    {
        return _imp->evaluate3(ir);
    }

    /*!
     * Construct the same observable as a regular observable
     */
    class TestRegularObservableProvider : public ParameterUser, public PrivateImplementationPattern<TestRegularObservableProvider>
    {
        public:
            TestRegularObservableProvider(const Parameters & parameters, const Options & options);
            ~TestRegularObservableProvider();

            // Observables
            double evaluate1(const double & q2) const;
            double evaluate2(const double & q2) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };

    const std::set<ReferenceName> TestRegularObservableProvider::references{};

    template <> struct Implementation<TestRegularObservableProvider>
    {
            UsedParameter m_B;

            double a;
            double b;

            Implementation(const Parameters & p, const Options & /* o */, ParameterUser & u) :
                m_B(p["mass::B_u"], u)
            {
                a = 2.0;
                b = m_B;
            }

            double
            evaluate1(const double & q2)
            {
                return b - a * q2;
            }

            double
            evaluate2(const double & q2)
            {
                return pow(q2, 2);
            }
    };

    TestRegularObservableProvider::TestRegularObservableProvider(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TestRegularObservableProvider>(new Implementation<TestRegularObservableProvider>(parameters, options, *this))
    {
    }

    TestRegularObservableProvider::~TestRegularObservableProvider() {}

    double
    TestRegularObservableProvider::evaluate1(const double & q2) const
    {
        return _imp->evaluate1(q2);
    }

    double
    TestRegularObservableProvider::evaluate2(const double & q2) const
    {
        return _imp->evaluate2(q2);
    }
} // namespace eos

class CacheableObservableTest : public TestCase
{
    public:
        CacheableObservableTest() :
            TestCase("cacheable_observable_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p   = Parameters::Defaults();
                p["mass::B_u"] = 5.27934;

                Options oo;

                TestCacheableObservableProvider tcop(p, oo);

                // Test the observable implementation
                TEST_CHECK_EQUAL(tcop.evaluate1(tcop.prepare(2.0)), 5.27934 - 2.0 * 2.0);
                TEST_CHECK_EQUAL(tcop.evaluate2(tcop.prepare(2.0)), 4.0);


                // Try to create a cacheable observable
                using TestCacheableObservable = class ConcreteCacheableObservable<TestCacheableObservableProvider, TestCacheableObservableProvider::IntermediateResult, double>;

                ObservablePtr cacheable_observable(new TestCacheableObservable("test::cacheable_observable1(q2)",
                                                                               p,
                                                                               Kinematics({
                                                                                   { "q2", 2.0 }
                }),
                                                                               Options(),
                                                                               &TestCacheableObservableProvider::prepare,
                                                                               &TestCacheableObservableProvider::evaluate1,
                                                                               std::make_tuple("q2")));

                TEST_CHECK_NEARLY_EQUAL(cacheable_observable->evaluate(), 5.27934 - 2.0 * 2.0, 1.e-5);

                // Try to add the observable to the cache...
                ObservableCache               cache(p);
                ObservableCache::ObservableId cacheable_observable_id;

                TEST_CHECK_NO_THROW(cacheable_observable_id = cache.add(cacheable_observable));

                // ..twice
                TEST_CHECK_EQUAL(cacheable_observable_id, cache.add(cacheable_observable));

                // Cache the observable by adding the same observable with a different name
                ObservablePtr                 cacheable_observable2(new TestCacheableObservable("test::cacheable_observable2(q2)",
                                                                                p,
                                                                                Kinematics({
                                                                                    { "q2", 2.0 }
                }),
                                                                                Options(),
                                                                                &TestCacheableObservableProvider::prepare,
                                                                                &TestCacheableObservableProvider::evaluate1,
                                                                                std::make_tuple("q2")));
                ObservableCache::ObservableId cacheable_observable2_id;
                unsigned                      cache_size = cache.size();

                TEST_CHECK_NO_THROW(cacheable_observable2_id = cache.add(cacheable_observable2));
                TEST_CHECK_EQUAL(cache.size(), cache_size + 1);


                // Create a regular observable
                using TestRegularObservable = class ConcreteObservable<TestRegularObservableProvider, double>;

                ObservablePtr                 regular_observable(new TestRegularObservable("test::regular_observable(q2)",
                                                                           p,
                                                                           Kinematics({
                                                                               { "q2", 2.0 }
                }),
                                                                           Options(),
                                                                           &TestRegularObservableProvider::evaluate1,
                                                                           std::make_tuple("q2")));
                ObservableCache::ObservableId regular_observable_id [[maybe_unused]];

                TEST_CHECK_NO_THROW(regular_observable_id = cache.add(regular_observable));


                // Add a third, different, cacheable observable
                ObservablePtr                 cacheable_observable3(new TestCacheableObservable("test::cacheable_observable3(q2)",
                                                                                p,
                                                                                Kinematics({
                                                                                    { "q2", 6.0 }
                }),
                                                                                Options(),
                                                                                &TestCacheableObservableProvider::prepare,
                                                                                &TestCacheableObservableProvider::evaluate2,
                                                                                std::make_tuple("q2")));
                ObservableCache::ObservableId cacheable_observable3_id [[maybe_unused]];

                TEST_CHECK_NO_THROW(cacheable_observable3_id = cache.add(cacheable_observable3));


                // Add a cacheable observable that shares its result type and its kinematics with the
                // first one, but obtains the intermediate result from a different prepare function
                ObservablePtr                 shifted_observable(new TestCacheableObservable("test::shifted_observable(q2)",
                                                                             p,
                                                                             Kinematics({
                                                                                 { "q2", 2.0 }
                }),
                                                                             Options(),
                                                                             &TestCacheableObservableProvider::prepare_shifted,
                                                                             &TestCacheableObservableProvider::evaluate1,
                                                                             std::make_tuple("q2")));
                ObservableCache::ObservableId shifted_observable_id;

                TEST_CHECK_NO_THROW(shifted_observable_id = cache.add(shifted_observable));


                // Add a cacheable observable that uses a second intermediate result type
                using TestAlternativeObservable = class ConcreteCacheableObservable<TestCacheableObservableProvider, TestCacheableObservableProvider::AlternativeResult, double>;

                ObservablePtr                 alternative_observable(new TestAlternativeObservable("test::alternative_observable(q2)",
                                                                                   p,
                                                                                   Kinematics({
                                                                                       { "q2", 2.0 }
                }),
                                                                                   Options(),
                                                                                   &TestCacheableObservableProvider::prepare_alternative,
                                                                                   &TestCacheableObservableProvider::evaluate3,
                                                                                   std::make_tuple("q2")));
                ObservableCache::ObservableId alternative_observable_id;

                TEST_CHECK_NO_THROW(alternative_observable_id = cache.add(alternative_observable));


                // Test cache evaluation
                TEST_CHECK_NO_THROW(cache.update());
                TEST_CHECK_EQUAL(cache[cacheable_observable_id], cache[cacheable_observable2_id]);
                TEST_CHECK_NEARLY_EQUAL(cache[shifted_observable_id], 5.27934 - 2.0 * 3.0, 1.e-5);
                TEST_CHECK_NEARLY_EQUAL(cache[alternative_observable_id], 5.27934 * 2.0, 1.e-5);

                // Test cache cloning
                ObservableCache cache2(p);
                TEST_CHECK_NO_THROW(cache2 = cache.clone(p));
            }
        }
} cacheable_observable_test;
