/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2024 Danny van Dyk
 *
 * Based upon 'paludis/util/log.hh', which is
 *
 *   Copyright (c) 2006-2011 Ciaran McCreesh
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

#ifndef EOS_GUARD_EOS_UTILS_LOG_HH
#define EOS_GUARD_EOS_UTILS_LOG_HH 1

#include <eos/utils/instantiation_policy.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/stringify.hh>

#include <functional>

namespace eos
{
    /*!
     * Categories for the severity of log messages
     */
    enum LogLevel
    {
        ll_silent,        ///< do not print any error message
        ll_error,         ///< only print error messages
        ll_warning,       ///< also print warning messages
        ll_success,       ///< also print success messages
        ll_completed,     ///< also print completion messages
        ll_inprogress,    ///< also print in-progress messages
        ll_informational, ///< also print informational messages
        ll_debug,         ///< also print debug messages
        ll_last
    };

    /*!
     * (De)stringification of LogLevel
     */
    ///@{
    std::ostream & operator<< (std::ostream & lhs, const LogLevel & rhs);
    std::istream & operator>> (std::istream & lhs, LogLevel & rhs);
    ///@}

    // forward declaration for use in Log
    class LogMessageHandler;

    /*!
     * Facility to emit messages to the user, with a filter according to interest and urgency
     * of the messages at hand.
     */
    class Log : public InstantiationPolicy<Log, Singleton>, public PrivateImplementationPattern<Log>
    {
        private:
            ///@name Basic Functions
            ///@{
            /// Constructor.
            Log();
            ///@}

            void _message(const std::string &, const LogLevel &, const std::string &);

        public:
            friend class LogMessageHandler;
            friend class InstantiationPolicy<Log, Singleton>;

            ///@name Basic Functions
            ///@{
            /// Destructor.
            ~Log();
            ///@}

            ///@name Access
            ///@{

            /// Get the current log level
            const LogLevel & get_log_level() const;

            /*!
             * Set the log level.
             */
            void set_log_level(const LogLevel &);

            /*!
             * Set the output stream.
             */
            void set_log_stream(std::ostream *);

            /*!
             * Set the program's name.
             */
            void set_program_name(const std::string &);

            /*!
             * Register a callback hook with the Log class.
             */
            void register_callback(const std::function<void(const std::string &, const LogLevel &, const std::string &)> &);

            /*!
             * Return a stream-like object to which message parts can be
             * appended via its overloaded operator<<. The message will be
             * completed upon destruction of the return value.
             */
            LogMessageHandler message(const std::string & id, const LogLevel & log_level) __attribute__((warn_unused_result));

            /*!
             * Class to be used as a RAII guard for a one-time log message.
             *
             * The message will be logged upon construction this object. Further messages can be avoided by
             * creating this object as a static const object in the scope of the message issuer.
             */
            class OneTimeMessage
            {
                public:
                    OneTimeMessage(const std::string & id, const LogLevel & log_level, const std::string & message);
            };
            friend class OneTimeMessage;
            ///@}
    };

    /*!
     * Proxy returned by Log::message.
     */
    class LogMessageHandler
    {
        private:
            Log *       _log;
            LogLevel    _log_level;
            std::string _id;
            std::string _message;

            ///@name Basic Functions
            ///@{
            /// Copy constructor.
            LogMessageHandler(const LogMessageHandler &);

            /// Constructor.
            LogMessageHandler(Log * const log, const LogLevel & log_level, const std::string & id);

            /// Assignment operator.
            void operator= (const LogMessageHandler &);
            ///@}

            void _append(const std::string & s);

        public:
            friend LogMessageHandler Log::message(const std::string &, const LogLevel &);

            ///@name Basic Function
            ///@{
            /// Destructor.
            ~LogMessageHandler();

            ///@}

            /*!
             * Append to our message.
             */
            template <typename T_>
            LogMessageHandler &
            operator<< (const T_ & t)
            {
                _append(stringify(t));

                return *this;
            }
    };
} // namespace eos

#endif
