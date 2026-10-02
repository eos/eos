/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/observable.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/observable_stub.hh>

#include <test/test.hh>

#include <set>
#include <stdexcept>

using namespace test;
using namespace eos;

namespace
{
    // A minimal Observable that always throws on evaluate(), used to exercise
    // ObservableCache::update()'s exception policy.
    class ThrowingObservable : public Observable
    {
        public:
            enum class Mode
            {
                eos_exception,
                std_exception
            };

        private:
            QualifiedName _name;
            Parameters    _parameters;
            Kinematics    _kinematics;
            Mode          _mode;

        public:
            ThrowingObservable(const Parameters & parameters, const QualifiedName & name, Mode mode) :
                _name(name),
                _parameters(parameters),
                _mode(mode)
            {
            }

            virtual ~ThrowingObservable() {}

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual double
            evaluate() const
            {
                if (Mode::eos_exception == _mode)
                {
                    throw InternalError("ThrowingObservable::evaluate(): deliberate failure");
                }

                throw std::runtime_error("ThrowingObservable::evaluate(): deliberate failure");
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Options
            options()
            {
                return Options();
            }

            virtual ObservablePtr
            clone() const
            {
                return ObservablePtr(new ThrowingObservable(_parameters.clone(), _name, _mode));
            }

            virtual ObservablePtr
            clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ThrowingObservable(parameters, _name, _mode));
            }
    };
} // namespace

class ObservableCacheEqualityTest : public TestCase
{
    public:
        ObservableCacheEqualityTest() :
            TestCase("observable_cache_equality_test")
        {
        }

        virtual void
        run() const
        {
            Parameters parameters = Parameters::Defaults();

            ObservableCache cache1(parameters);
            ObservableCache cache2 = cache1;
            ObservableCache cache3(parameters);

            TEST_CHECK(cache1 == cache2);
            TEST_CHECK(! (cache1 == cache3));
        }
} observable_cache_equality_test;

class ObservableCacheGenerationTest : public TestCase
{
    public:
        ObservableCacheGenerationTest() :
            TestCase("observable_cache_generation_test")
        {
        }

        virtual void
        run() const
        {
            Parameters      parameters = Parameters::Defaults();
            ObservableCache cache(parameters);

            cache.add(ObservablePtr(new ObservableStub(parameters, "mass::c")));

            const ObservableCache::Generation initial = cache.generation();
            TEST_CHECK(ObservableCache::Generation() != initial);
            TEST_CHECK(cache.generation() == initial);

            cache.update();
            const ObservableCache::Generation first = cache.generation();
            TEST_CHECK(first != initial);
            TEST_CHECK(cache.generation() == first);

            cache.update();
            TEST_CHECK(cache.generation() != first);

            // distinct caches never share a generation, not even before their first update()
            ObservableCache other(parameters);
            TEST_CHECK(other.generation() != initial);
        }
} observable_cache_generation_test;

class ObservableCacheExceptionPolicyTest : public TestCase
{
    public:
        ObservableCacheExceptionPolicyTest() :
            TestCase("observable_cache_exception_policy_test")
        {
        }

        virtual void
        run() const
        {
            // an eos::Exception thrown by a regular observable surfaces on the calling thread,
            // and the generation counter does not advance
            {
                Parameters      parameters = Parameters::Defaults();
                ObservableCache cache(parameters);

                cache.add(ObservablePtr(new ThrowingObservable(parameters, "mass::c", ThrowingObservable::Mode::eos_exception)));

                const ObservableCache::Generation initial = cache.generation();
                TEST_CHECK_THROWS(eos::Exception, cache.update());
                TEST_CHECK(cache.generation() == initial);
            }

            // a non-eos::Exception (e.g. from a GSL/Boost wrapper or an external observable) is
            // caught rather than escaping into the thread pool and terminating the process
            {
                Parameters      parameters = Parameters::Defaults();
                ObservableCache cache(parameters);

                cache.add(ObservablePtr(new ThrowingObservable(parameters, "mass::c", ThrowingObservable::Mode::std_exception)));

                const ObservableCache::Generation initial = cache.generation();
                TEST_CHECK_THROWS(std::runtime_error, cache.update());
                TEST_CHECK(cache.generation() == initial);
            }

            // with two failing observables in one batch, update() throws exactly once and the
            // process survives; the cache remains usable afterwards
            {
                Parameters      parameters = Parameters::Defaults();
                ObservableCache cache(parameters);

                std::vector<ObservablePtr> batch;
                batch.push_back(ObservablePtr(new ThrowingObservable(parameters, "mass::c", ThrowingObservable::Mode::eos_exception)));
                batch.push_back(ObservablePtr(new ThrowingObservable(parameters, "mass::b(MSbar)", ThrowingObservable::Mode::std_exception)));
                cache.add_batch(std::move(batch));

                const ObservableCache::Generation initial = cache.generation();
                TEST_CHECK_THROWS(std::exception, cache.update());
                TEST_CHECK(cache.generation() == initial);

                // the cache survives and a fresh update() attempt behaves the same way
                TEST_CHECK_THROWS(std::exception, cache.update());
                TEST_CHECK(cache.generation() == initial);
            }
        }
} observable_cache_exception_policy_test;

class ObservableCacheCloneTest : public TestCase
{
    public:
        ObservableCacheCloneTest() :
            TestCase("observable_cache_clone_test")
        {
        }

        virtual void
        run() const
        {
            Parameters      parameters = Parameters::Defaults();
            ObservableCache cache(parameters);

            const auto id = cache.add(ObservablePtr(new ObservableStub(parameters, "mass::c")));

            std::vector<ObservablePtr> batch;
            batch.push_back(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")));
            const auto batch_id = cache.add_batch(std::move(batch));
            cache.update();

            // ObservableIds remain valid for the clone, while batches are left to their owners
            Parameters      clone_parameters = parameters.clone();
            ObservableCache clone            = cache.clone(clone_parameters);
            TEST_CHECK_EQUAL(clone[id], cache[id]);
            TEST_CHECK_THROWS(InternalError, clone[batch_id]);
            TEST_CHECK_EQUAL(cache[batch_id].size(), 1u);
        }
} observable_cache_clone_test;

class ObservableCacheUsedParametersTest : public TestCase
{
    public:
        ObservableCacheUsedParametersTest() :
            TestCase("observable_cache_used_parameters_test")
        {
        }

        virtual void
        run() const
        {
            Parameters      parameters = Parameters::Defaults();
            ObservableCache cache(parameters);

            // register through a second handle, as a LogLikelihoodBlock does; both share one implementation
            {
                ObservableCache handle = cache;
                handle.add(ObservablePtr(new ObservableStub(parameters, "mass::c")));

                std::vector<ObservablePtr> batch;
                batch.push_back(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")));
                handle.add_batch(std::move(batch));
            }

            const std::set<Parameter::Id> used(cache.begin_used_parameter_ids(), cache.end_used_parameter_ids());
            const std::set<Parameter::Id> expected{ parameters["mass::c"].id(), parameters["mass::b(MSbar)"].id() };
            TEST_CHECK(used == expected);
        }
} observable_cache_used_parameters_test;
