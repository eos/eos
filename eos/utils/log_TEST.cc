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
#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>

using namespace test;
using namespace eos;

class LogLevelTest :
    public TestCase
{
    public:
        LogLevelTest() :
            TestCase("log_level_test")
        {
        }

        virtual void run() const
        {
            // stringification
            {
                TEST_CHECK_EQUAL("silent",        stringify(ll_silent));
                TEST_CHECK_EQUAL("error",         stringify(ll_error));
                TEST_CHECK_EQUAL("warning",       stringify(ll_warning));
                TEST_CHECK_EQUAL("informational", stringify(ll_informational));
                TEST_CHECK_EQUAL("debug",         stringify(ll_debug));
                TEST_CHECK_THROWS(InternalError,  stringify(ll_last));
            }

            // destringification
            {
                // silent
                {
                    std::stringstream ss("silent");
                    LogLevel log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_silent);
                }

                // error
                {
                    std::stringstream ss("error");
                    LogLevel log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_error);
                }

                // warning
                {
                    std::stringstream ss("warning");
                    LogLevel log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_warning);
                }

                // informational
                {
                    std::stringstream ss("informational");
                    LogLevel log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_informational);
                }

                // debug
                {
                    std::stringstream ss("debug");
                    LogLevel log_level;

                    ss >> log_level;

                    TEST_CHECK_EQUAL(log_level, ll_debug);
                }

                // last
                {
                    std::stringstream ss("last");
                    LogLevel log_level;

                    TEST_CHECK_THROWS(InternalError, ss >> log_level);
                }
            }
        }
} log_level_test;
