/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
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

#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>

#include <test/test.hh>

#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

using namespace test;
using namespace eos;

class LogLevelTest : public TestCase
{
    public:
        LogLevelTest() :
            TestCase("log_level_test")
        {
        }

        virtual void
        run() const
        {
            // stringification
            {
                TEST_CHECK_EQUAL("silent", stringify(ll_silent));
                TEST_CHECK_EQUAL("error", stringify(ll_error));
                TEST_CHECK_EQUAL("warning", stringify(ll_warning));
                TEST_CHECK_EQUAL("success", stringify(ll_success));
                TEST_CHECK_EQUAL("completed", stringify(ll_completed));
                TEST_CHECK_EQUAL("inprogress", stringify(ll_inprogress));
                TEST_CHECK_EQUAL("informational", stringify(ll_informational));
                TEST_CHECK_EQUAL("debug", stringify(ll_debug));
                TEST_CHECK_THROWS(InternalError, stringify(ll_last));
            }

            // destringification
            {
                // silent
                {
                    std::stringstream ss("silent");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_silent);
                }

                // error
                {
                    std::stringstream ss("error");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_error);
                }

                // warning
                {
                    std::stringstream ss("warning");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_warning);
                }

                // success
                {
                    std::stringstream ss("success");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_success);
                }

                // completed
                {
                    std::stringstream ss("completed");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_completed);
                }

                // inprogress
                {
                    std::stringstream ss("inprogress");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_inprogress);
                }

                // informational
                {
                    std::stringstream ss("informational");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_informational);
                }

                // debug
                {
                    std::stringstream ss("debug");
                    LogLevel          log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_debug);
                }

                // last
                {
                    std::stringstream ss("last");
                    LogLevel          log_level;

                    TEST_CHECK_THROWS(InternalError, ss >> log_level);
                }
            }
        }
} log_level_test;

class LogOneTimeMessageTest : public TestCase
{
    public:
        LogOneTimeMessageTest() :
            TestCase("log_one_time_message_test")
        {
        }

        // the callback outlives run() and may be invoked concurrently, so it must neither capture
        // anything local to run() nor record a message without holding our mutex
        mutable std::mutex                                                  mutex;
        mutable std::vector<std::tuple<std::string, LogLevel, std::string>> messages;

        unsigned
        number_of_messages() const
        {
            std::unique_lock<std::mutex> l(mutex);

            return messages.size();
        }

        virtual void
        run() const
        {
            static const std::string id_of_interest = "test-one-time-message";

            // register callback
            std::function<void(const std::string &, const LogLevel &, const std::string &)> callback =
                    [this](const std::string & id, const LogLevel & level, const std::string & message)
            {
                if (id_of_interest != id)
                {
                    return;
                }

                std::unique_lock<std::mutex> l(mutex);

                messages.push_back(std::make_tuple(id, level, message));
            };
            Log::instance()->register_callback(callback);

            // first message
            {
                // no messages yet
                TEST_CHECK_EQUAL(0u, number_of_messages());

                // emit message
                static const Log::OneTimeMessage first_emission(id_of_interest, ll_informational, "This is a test message.");

                // one message
                TEST_CHECK_EQUAL(1u, number_of_messages());

                // check message
                TEST_CHECK_EQUAL(id_of_interest, std::get<0>(messages.back()));
                TEST_CHECK_EQUAL(ll_informational, std::get<1>(messages.back()));
                TEST_CHECK_EQUAL("This is a test message. (Further messages of this type will be suppressed.)", std::get<2>(messages.back()));

                // try emitting another message with the same id
                static const Log::OneTimeMessage second_emission(id_of_interest, ll_informational, "This is a test message.");

                // still one message
                TEST_CHECK_EQUAL(1u, number_of_messages());
            }
        }
} one_time_message_test;

class LogCallbackLockingTest : public TestCase
{
    public:
        LogCallbackLockingTest() :
            TestCase("log_callback_locking_test")
        {
        }

        // the callback outlives run(), so it must not capture anything local to it
        mutable std::mutex              mutex;
        mutable std::condition_variable condition;
        mutable bool                    within_callback      = false;
        mutable bool                    other_emitted        = false;
        mutable bool                    emitted_concurrently = false;

        static const std::string &
        id()
        {
            static const std::string result = "test-callback-locking";

            return result;
        }

        virtual void
        run() const
        {
            static const auto timeout = std::chrono::seconds(5);

            // this callback holds up the emitting thread until a second thread has emitted a
            // message of its own; that second message needs the Log's lock, which must therefore
            // not be held while we are here
            std::function<void(const std::string &, const LogLevel &, const std::string &)> callback = [this](const std::string & i, const LogLevel &, const std::string &)
            {
                if (id() != i)
                {
                    return;
                }

                std::unique_lock<std::mutex> l(mutex);

                within_callback = true;
                condition.notify_all();

                // the second thread reaches us only if we do not hold the Log's lock; waiting in
                // vain means that it is blocked on that lock until we return
                emitted_concurrently = condition.wait_for(l, timeout, [this]() { return other_emitted; });
            };
            Log::instance()->register_callback(callback);

            std::thread other([this]()
            {
                {
                    std::unique_lock<std::mutex> l(mutex);

                    if (! condition.wait_for(l, timeout, [this]() { return within_callback; }))
                    {
                        return;
                    }
                }

                Log::instance()->message(id() + "-other", ll_debug) << "emitted while the first thread is within the callback";

                std::unique_lock<std::mutex> l(mutex);

                other_emitted = true;
                condition.notify_all();
            });

            Log::instance()->message(id(), ll_debug) << "emitted from the first thread";

            other.join();

            TEST_CHECK(other_emitted);
            TEST_CHECK(emitted_concurrently);
        }
} log_callback_locking_test;
