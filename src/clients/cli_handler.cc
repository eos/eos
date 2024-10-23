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

#include "cli_handler.hh"

#include <eos/utils/indirect-iterator-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include "cli_dumper.hh"
#include "cli_error.hh"
#include "cli_visitor.hh"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>

namespace eos
{
    /**
     * Imp data for Handler.
     */
    template <> struct Implementation<cli::Handler>
    {
            std::list<cli::Section *>                      sections;
            std::list<std::string>                         parameters;
            std::list<std::string>                         usage_lines;
            std::list<std::pair<std::string, std::string>> environment_lines;
            std::list<std::pair<std::string, std::string>> example_lines;
            std::list<std::string>                         notes;
            std::list<std::string>                         descriptions;
            std::list<std::pair<std::string, int>>         see_alsos;

            std::map<std::string, cli::Option *> longopts;
            std::map<char, cli::Option *>        shortopts;

            std::shared_ptr<cli::Section> main_options_section;

            Implementation() {}
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::SectionsConstIteratorTag>
    {
            using UnderlyingIterator = IndirectIterator<std::list<cli::Section *>::const_iterator>;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::ParametersConstIteratorTag>
    {
            using UnderlyingIterator = std::list<std::string>::const_iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::NotesIteratorTag>
    {
            using UnderlyingIterator = std::list<std::string>::const_iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::DescriptionLineConstIteratorTag>
    {
            using UnderlyingIterator = std::list<std::string>::const_iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::ExamplesConstIteratorTag>
    {
            using UnderlyingIterator = std::list<std::pair<std::string, std::string>>::const_iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::UsageLineConstIteratorTag>
    {
            using UnderlyingIterator = std::list<std::string>::const_iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::ArgsIteratorTag>
    {
            using UnderlyingIterator = std::list<std::string>::iterator;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Handler::SeeAlsoConstIteratorTag>
    {
            using UnderlyingIterator = std::list<std::pair<std::string, int>>::const_iterator;
    };

    namespace cli
    {
        Handler::Handler() :
            PrivateImplementationPattern<cli::Handler>(new Implementation<cli::Handler>())
        {
        }

        Handler::~Handler() = default;

        void
        Handler::add_usage_line(const std::string & l)
        {
            _imp->usage_lines.push_back(l);
        }

        void
        Handler::add_example(const std::string & e, const std::string & desc)
        {
            _imp->example_lines.push_back(std::make_pair(e, desc));
        }

        void
        Handler::add_note(const std::string & e)
        {
            _imp->notes.push_back(e);
        }

        void
        Handler::add_see_also(const std::string & e, const int s)
        {
            _imp->see_alsos.push_back(std::make_pair(e, s));
        }

        void
        Handler::add(Section * const s)
        {
            _imp->sections.push_back(s);
        }

        void
        Handler::run(const int argc, const char * const * const argv, const std::string & /*client*/)
        {
            std::list<std::string> args;
            std::copy(&argv[1], &argv[argc], std::back_inserter(args));
            ArgsIterator argit(args.begin()), arge(args.end());

            for (; argit != arge; ++argit)
            {
                std::string arg = *argit;

                if (arg == "--")
                {
                    throw BadArgument("--");
                }
                else if (0 == arg.compare(0, 2, "--"))
                {
                    arg.erase(0, 2);
                    std::map<std::string, cli::Option *>::iterator it = _imp->longopts.find(arg);
                    if (it == _imp->longopts.end())
                    {
                        if (0 != arg.compare(0, 3, "no-"))
                        {
                            throw BadArgument("--" + arg);
                        }
                        arg.erase(0, 3);
                        it = _imp->longopts.find(arg);
                        if (it == _imp->longopts.end())
                        {
                            throw BadArgument("--no-" + arg);
                        }
                        if (! it->second->can_be_negated())
                        {
                            throw BadArgument("--no-" + arg);
                        }

                        std::string remaining_chars;
                        Visitor     visitor(&argit, arge, remaining_chars, true);
                        it->second->accept(visitor);
                    }
                    else
                    {
                        std::string remaining_chars;
                        Visitor     visitor(&argit, arge, remaining_chars, false);
                        it->second->accept(visitor);
                    }
                }
                else if (arg[0] == '-' || arg[0] == '+')
                {
                    bool negate(arg[0] == '+');
                    arg.erase(0, 1);
                    for (std::string::iterator c = arg.begin(); c != arg.end(); ++c)
                    {
                        bool        maybe_second_char_used(false);
                        std::string remaining_chars;
                        if (arg.length() >= 2 && c == arg.begin())
                        {
                            maybe_second_char_used = true;
                            remaining_chars        = arg.substr(1);
                        }

                        std::map<char, Option *>::iterator it = _imp->shortopts.find(*c);
                        if (it == _imp->shortopts.end())
                        {
                            throw BadArgument(std::string(negate ? "+" : "-") + *c);
                        }
                        if (negate && ! it->second->can_be_negated())
                        {
                            throw BadArgument(std::string("+") + *c);
                        }

                        Visitor visitor(&argit, arge, remaining_chars, negate);
                        it->second->accept(visitor);

                        if (maybe_second_char_used && remaining_chars.empty())
                        {
                            break;
                        }
                    }
                }
                else
                {
                    _imp->parameters.push_back(arg);
                }
            }

            _imp->parameters.insert(_imp->parameters.end(), argit, ArgsIterator(args.end()));

            post_run();
        }

        void
        Handler::post_run()
        {
        }

        void
        Handler::dump_to_stream(std::ostream & s) const
        {
            Dumper dump(s);
            for (SectionsConstIterator a(begin_args_sections()), a_end(end_args_sections()); a != a_end; ++a)
            {
                for (Section::GroupsConstIterator g(a->begin()), g_end(a->end()); g != g_end; ++g)
                {
                    s << g->name() << ":" << std::endl;

                    std::for_each(indirect_iterator(g->begin()), indirect_iterator(g->end()), accept_visitor(dump));

                    s << std::endl;
                }
            }
        }

        Handler::ParametersConstIterator
        Handler::begin_parameters() const
        {
            return ParametersConstIterator(_imp->parameters.begin());
        }

        Handler::ParametersConstIterator
        Handler::end_parameters() const
        {
            return ParametersConstIterator(_imp->parameters.end());
        }

        bool
        Handler::empty() const
        {
            return _imp->parameters.empty();
        }

        std::ostream &
        operator<< (std::ostream & s, const Handler & h)
        {
            h.dump_to_stream(s);
            return s;
        }

        void
        Handler::add_option(Option * const opt, const std::string & long_name, const char short_name)
        {
            if (! _imp->longopts.insert(std::make_pair(long_name, opt)).second)
            {
                throw InternalError("duplicate long name '" + long_name + "'");
            }

            _imp->longopts[long_name] = opt;
            if (short_name != '\0')
            {
                if (! _imp->shortopts.insert(std::make_pair(short_name, opt)).second)
                {
                    throw InternalError("duplicate short name '" + stringify(short_name) + "'");
                }
            }
        }

        void
        Handler::remove_option(const std::string & long_name, const char short_name)
        {
            _imp->longopts.erase(long_name);
            if (short_name != '\0')
            {
                _imp->shortopts.erase(short_name);
            }
        }

        Handler::UsageLineConstIterator
        Handler::begin_usage_lines() const
        {
            return UsageLineConstIterator(_imp->usage_lines.begin());
        }

        Handler::UsageLineConstIterator
        Handler::end_usage_lines() const
        {
            return UsageLineConstIterator(_imp->usage_lines.end());
        }

        Handler::ExamplesConstIterator
        Handler::begin_examples() const
        {
            return ExamplesConstIterator(_imp->example_lines.begin());
        }

        Handler::ExamplesConstIterator
        Handler::end_examples() const
        {
            return ExamplesConstIterator(_imp->example_lines.end());
        }

        Handler::NotesIterator
        Handler::begin_notes() const
        {
            return NotesIterator(_imp->notes.begin());
        }

        Handler::NotesIterator
        Handler::end_notes() const
        {
            return NotesIterator(_imp->notes.end());
        }

        Handler::SeeAlsoConstIterator
        Handler::begin_see_alsos() const
        {
            return SeeAlsoConstIterator(_imp->see_alsos.begin());
        }

        Handler::SeeAlsoConstIterator
        Handler::end_see_alsos() const
        {
            return SeeAlsoConstIterator(_imp->see_alsos.end());
        }

        Handler::SectionsConstIterator
        Handler::begin_args_sections() const
        {
            return SectionsConstIterator(indirect_iterator(_imp->sections.cbegin()));
        }

        Handler::SectionsConstIterator
        Handler::end_args_sections() const
        {
            return SectionsConstIterator(indirect_iterator(_imp->sections.cend()));
        }

        Section *
        Handler::main_options_section()
        {
            if (! _imp->main_options_section)
            {
                _imp->main_options_section = std::make_shared<Section>(this, "Options");
            }
            return _imp->main_options_section.get();
        }

        Handler::DescriptionLineConstIterator
        Handler::begin_description_lines() const
        {
            return DescriptionLineConstIterator(_imp->descriptions.begin());
        }

        Handler::DescriptionLineConstIterator
        Handler::end_description_lines() const
        {
            return DescriptionLineConstIterator(_imp->descriptions.end());
        }

        void
        Handler::add_description_line(const std::string & e)
        {
            _imp->descriptions.push_back(e);
        }
    } // namespace cli

    template class WrappedForwardIterator<cli::Handler::ParametersConstIteratorTag, const std::string>;
    template class WrappedForwardIterator<cli::Handler::UsageLineConstIteratorTag, const std::string>;
    template class WrappedForwardIterator<cli::Handler::ExamplesConstIteratorTag, const std::pair<std::string, std::string>>;
    template class WrappedForwardIterator<cli::Handler::SectionsConstIteratorTag, const cli::Section>;
    template class WrappedForwardIterator<cli::Handler::NotesIteratorTag, const std::string>;
    template class WrappedForwardIterator<cli::Handler::DescriptionLineConstIteratorTag, const std::string>;
    template class WrappedForwardIterator<cli::Handler::ArgsIteratorTag, std::string>;
    template class WrappedForwardIterator<cli::Handler::SeeAlsoConstIteratorTag, const std::pair<std::string, int>>;

    namespace cli
    {
        DefaultHandler::DefaultHandler() :
            g_universal_options(main_options_section(), "Universal Options", "Universal options, common to all command-line clients."),
            a_help(&g_universal_options, "help", 'h', "display help message", false),
            a_log_level(&g_universal_options, "log-level", 'L'),
            a_version(&g_universal_options, "version", 'v', "display version information", false)
        {
            add_usage_line("[ --log-level level ]");
            add_usage_line("help [ --all ]");
        }
    } // namespace cli
} // namespace eos
