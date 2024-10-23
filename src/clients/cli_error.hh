/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_ERROR_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_ERROR_HH 1

#include <eos/utils/exception.hh>

#include <string>

namespace eos
{
    namespace cli
    {
        /**
         * Thrown if an invalid command line argument is provided.
         */
        class Error : public eos::Exception
        {
            protected:
                /**
                 * Constructor.
                 */
                Error(const std::string & message) noexcept;
        };

        /**
         * Thrown if an unrecognised command line argument is specified.
         */
        class BadArgument : public cli::Error
        {
            public:
                /**
                 * Constructor.
                 */
                BadArgument(const std::string & option) noexcept;
        };

        /**
         * Thrown if an invalid parameter is passed to a valid command line argument.
         */
        class BadValue : public cli::Error
        {
            public:
                /**
                 * Constructor
                 */
                BadValue(const std::string & option, const std::string & value) noexcept;
        };

        /**
         * Thrown if an argument is specified that needs a parameter,
         * but no parameter is given.
         */
        class MissingValue : public cli::Error
        {
            public:
                /**
                 * Constructor.
                 */
                MissingValue(const std::string & arg) noexcept;
        };

        /**
         * Thrown to signal that the help message needs to be displayed.
         */
        struct DoHelp
        {
                const std::string message;

                DoHelp(const std::string & m = "");
        };
    } // namespace cli
} // namespace eos

#endif
