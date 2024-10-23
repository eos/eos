/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2005-2011 Ciaran McCreesh
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

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_HANDLER_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_HANDLER_HH 1

#include <eos/utils/private_implementation_pattern.hh>

#include "cli_group.hh"
#include "cli_section.hh"
#include <iosfwd>
#include <memory>
#include <string>

namespace eos
{
    namespace cli
    {
        class Handler : public PrivateImplementationPattern<cli::Handler>
        {
                friend class cli::Section;
                friend std::ostream & operator<< (std::ostream &, const Handler &);

            protected:
                /**
                 * Add a new usage line.
                 */
                void add_usage_line(const std::string & l);

                /**
                 * Add a new example.
                 */
                void add_example(const std::string & e, const std::string & desc);

                /**
                 * Add a new note.
                 */
                void add_note(const std::string &);

                /**
                 * Add a new description.
                 */
                void add_description_line(const std::string & l);

                /**
                 * Add a 'see also' item.
                 */
                void add_see_also(const std::string &, int section);

                /**
                 * Add a new Section (called by the Section constructor).
                 */
                void add(Section * const);

                /**
                 * Dump, for --help output (called by operator<<).
                 */
                void dump_to_stream(std::ostream & s) const;

                /**
                 * Called after run(), for convenience. Does nothing.
                 */
                virtual void post_run();

            public:
                ///@name Basic operations
                ///@{

                Handler();

                virtual ~Handler();

                Handler(const Handler &) = delete;

                Handler & operator= (const Handler &) = delete;

                ///@}

                ///@name Iterate over our parameters (non - and -- switches and their values)
                ///@{

                struct ParametersConstIteratorTag;
                using ParametersConstIterator = WrappedForwardIterator<ParametersConstIteratorTag, const std::string>;

                ParametersConstIterator begin_parameters() const;
                ParametersConstIterator end_parameters() const;
#if 0
                IteratorRange<ParametersConstIterator> parameters() const {
                    return make_range(begin_parameters(), end_parameters());
                }
#endif
                bool empty() const;

                /**
                 * Add an Option instance.
                 */
                void add_option(cli::Option * const, const std::string & long_name, const char short_name = '\0');

                /**
                 * Remove an Option instance.
                 */
                void remove_option(const std::string & long_name, const char short_name = '\0');

                ///@name About our application (for documentation)
                ///@{

                /**
                 * What is our application name?
                 */
                virtual std::string app_name() const = 0;

                /**
                 * What is our application's Unix manual section?
                 */
                virtual std::string
                man_section() const
                {
                    return "1";
                }

                /**
                 * One line synopsis of what our application is.
                 */
                virtual std::string app_synopsis() const = 0;

                /**
                 * Long description of what our application is.
                 */
                virtual std::string app_description() const = 0;

                ///@}

                ///@name Iterate over our usage lines (for documentation)
                ///@{

                struct UsageLineConstIteratorTag;
                using UsageLineConstIterator = WrappedForwardIterator<UsageLineConstIteratorTag, const std::string>;

                UsageLineConstIterator begin_usage_lines() const;
                UsageLineConstIterator end_usage_lines() const;
#if 0
                IteratorRange<UsageLineConstIterator> usage_lines() const {
                    return make_range(begin_usage_lines(), end_usage_lines());
                }
#endif
                ///@}

                ///@name Iterate over our examples (for documentation)
                ///@{

                struct ExamplesConstIteratorTag;
                using ExamplesConstIterator = WrappedForwardIterator<ExamplesConstIteratorTag, const std::pair<std::string, std::string>>;

                ExamplesConstIterator begin_examples() const;
                ExamplesConstIterator end_examples() const;
#if 0
                IteratorRange<ExamplesConstIterator> examples() const {
                    return make_range(begin_examples(), end_examples());
                }
#endif
                ///@}

                ///@name Iterate over our sections
                ///@{

                struct SectionsConstIteratorTag;
                using SectionsConstIterator = WrappedForwardIterator<SectionsConstIteratorTag, const Section>;

                SectionsConstIterator begin_args_sections() const;
                SectionsConstIterator end_args_sections() const;
#if 0
                IteratorRange<SectionsConstIterator> args_sections() const {
                    return make_range(begin_args_sections(), end_args_sections());
                }
#endif
                /**
                 * The 'Options' section.
                 *
                 * Created if it does not exist.
                 */
                Section * main_options_section();

                ///@}

                ///@name Iterate over our notes
                ///@{

                struct NotesIteratorTag;
                using NotesIterator = WrappedForwardIterator<NotesIteratorTag, const std::string>;

                NotesIterator begin_notes() const;
                NotesIterator end_notes() const;
#if 0
                IteratorRange<NotesIterator> notes() const {
                    return make_range(begin_notes(), end_notes());
                }
#endif
                ///@}

                ///@name Iterate over our extra description lines (for documentation)
                ///@{

                struct DescriptionLineConstIteratorTag;
                using DescriptionLineConstIterator = WrappedForwardIterator<DescriptionLineConstIteratorTag, const std::string>;

                DescriptionLineConstIterator begin_description_lines() const;
                DescriptionLineConstIterator end_description_lines() const;
#if 0
                IteratorRange<DescriptionLineConstIterator> description_lines() const {
                    return make_range(begin_description_lines(), end_description_lines());
                }
#endif
                ///@}

                ///@name Iterate over our 'see also' lines
                ///@since 0.48.2
                ///@{

                struct SeeAlsoConstIteratorTag;
                using SeeAlsoConstIterator = WrappedForwardIterator<SeeAlsoConstIteratorTag, const std::pair<std::string, int>>;

                SeeAlsoConstIterator begin_see_alsos() const;
                SeeAlsoConstIterator end_see_alsos() const;
#if 0
                IteratorRange<SeeAlsoConstIterator> see_alsos() const {
                    return make_range(begin_see_alsos(), end_see_alsos());
                }
#endif
                ///@}

                ///@name For use by ArgsVisitor
                ///@{

                struct ArgsIteratorTag;
                using ArgsIterator = WrappedForwardIterator<ArgsIteratorTag, std::string>;

                ///@}

                /**
                 * Parse command line arguments.
                 */
                void run(const int argc, const char * const * const argv, const std::string & client);
        };

        /**
         * Output an Handler to an ostream, for --help output.
         */
        std::ostream & operator<< (std::ostream &, const Handler &);
    } // namespace cli

    extern template class WrappedForwardIterator<cli::Handler::ParametersConstIteratorTag, const std::string>;
    extern template class WrappedForwardIterator<cli::Handler::UsageLineConstIteratorTag, const std::string>;
    extern template class WrappedForwardIterator<cli::Handler::ExamplesConstIteratorTag, const std::pair<std::string, std::string>>;
    extern template class WrappedForwardIterator<cli::Handler::SectionsConstIteratorTag, const cli::Section>;
    extern template class WrappedForwardIterator<cli::Handler::DescriptionLineConstIteratorTag, const std::string>;
    extern template class WrappedForwardIterator<cli::Handler::NotesIteratorTag, const std::string>;
    extern template class WrappedForwardIterator<cli::Handler::ArgsIteratorTag, std::string>;
    extern template class WrappedForwardIterator<cli::Handler::SeeAlsoConstIteratorTag, const std::pair<std::string, int>>;

    namespace cli
    {
        /**
         * The default command line Handler.
         *
         * Knows about --help, --log-level, and --version.
         */
        class DefaultHandler : public Handler
        {
            public:
                // universal options
                cli::Group       g_universal_options;
                cli::SwitchArg   a_help;
                cli::LogLevelArg a_log_level;
                cli::SwitchArg   a_version;

                DefaultHandler();

                virtual ~DefaultHandler() = default;
        };
    } // namespace cli
} // namespace eos

#endif
