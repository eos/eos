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

#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/instantiation_policy.hh>

#include <test/test.hh>

#include <atomic>
#include <thread>
#include <vector>

using namespace test;
using namespace eos;

namespace
{
    std::atomic<unsigned> constructions(0);

    class Counted : public InstantiationPolicy<Counted, Singleton>
    {
        public:
            unsigned payload;

            Counted() :
                payload(4711u)
            {
                constructions.fetch_add(1u);
            }
    };
} // namespace

class InstantiationPolicySingletonTest : public TestCase
{
    public:
        InstantiationPolicySingletonTest() :
            TestCase("instantiation_policy_singleton_test")
        {
        }

        virtual void
        run() const
        {
            static const unsigned number_of_threads = 16;

            // every thread must observe the one instance, fully constructed
            std::vector<std::thread> threads;
            std::vector<Counted *>   instances(number_of_threads, nullptr);
            std::atomic<unsigned>    ready(0);

            for (unsigned i = 0; i < number_of_threads; ++i)
            {
                threads.emplace_back([&instances, &ready, i]()
                {
                    // line the threads up, so that they race for the first call
                    ready.fetch_add(1u);
                    while (ready.load() < number_of_threads)
                    {
                    }

                    instances[i] = Counted::instance();
                });
            }

            for (auto & thread : threads)
            {
                thread.join();
            }

            TEST_CHECK_EQUAL(1u, constructions.load());

            for (unsigned i = 0; i < number_of_threads; ++i)
            {
                TEST_CHECK(nullptr != instances[i]);
                TEST_CHECK_EQUAL(instances[0], instances[i]);
                TEST_CHECK_EQUAL(4711u, instances[i]->payload);
            }
        }
} instantiation_policy_singleton_test;
