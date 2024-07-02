/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2023 Danny van Dyk
 * Copyright (c) 2016 Rafael Silva Coutinho
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

#include <eos/utils/qualified-name.hh>

namespace eos
{
    namespace qnp
    {
        Prefix::Prefix(const std::string & prefix) :
            _prefix(prefix)
        {
            if (prefix.empty())
            {
                throw QualifiedNameSyntaxError("A qualified name's prefix part must not be empty");
            }

            // PREFIX := ['a'-'z', 'A'-'Z', '0'-'9', '<', '>', '^', '_', '*', '+', '-', '(', ')']
            static const char * valid_prefix_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789"
                    "<>^_*+-()";

            auto pos = prefix.find_first_not_of(valid_prefix_characters);

            if (std::string::npos != pos)
            {
                throw QualifiedNameSyntaxError("'" + prefix + "' is not a valid prefix part: Character '" + prefix[pos] + "' may not be used");
            }

            auto pos_opened = prefix.find("<<");
            auto pos_closed = prefix.find(">>");

            if ((std::string::npos != pos_opened) || (std::string::npos != pos_closed))
            {
                throw QualifiedNameSyntaxError("'" + prefix + "' is not a valid prefix part: Neither '<<' nor '>>' may be used");
            }
        }

        Name::Name(const std::string & name) :
            _name(name)
        {
            if (name.empty())
            {
                throw QualifiedNameSyntaxError("A qualified name's name part must not be empty");
            }

            // NAME := ['a'-'z', 'A'-'Z', '0'-'9', '(', ')', '[', ']', '{', '}', '|'
            //          '\'', ',', '.', '/', '^', '_', '*', '+', '-', '=']
            static const char * valid_name_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789"
                    "()[]{}|',./^_*+-=";

            auto pos = name.find_first_not_of(valid_name_characters);

            if (std::string::npos != pos)
            {
                throw QualifiedNameSyntaxError("'" + name + "' is not a valid name part: Character '" + name[pos] + "' may not be used");
            }

            auto pos_opened = name.find("[[");
            auto pos_closed = name.find("]]");

            if ((std::string::npos != pos_opened) || (std::string::npos != pos_closed))
            {
                throw QualifiedNameSyntaxError("'" + name + "' is not a valid prefix part: Neither '[[' nor ']]' may be used");
            }
        }

        Suffix::Suffix() :
            _suffix()
        {
        }

        Suffix::Suffix(const std::string & suffix) :
            _suffix(suffix)
        {
            // SUFFIX := ['a'-'z', 'A'-'Z', '0'-'9', '.', ':', '-', '(', ')']
            static const char * valid_suffix_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789"
                    ".:-+()";

            auto pos = suffix.find_first_not_of(valid_suffix_characters);

            if (std::string::npos != pos)
            {
                throw QualifiedNameSyntaxError("'" + suffix + "' is not a valid suffix part: Character '" + suffix[pos] + "' may not be used");
            }
        }

        OptionKey::OptionKey(const std::string & key) :
            _key(key)
        {
            if (key.empty())
            {
                throw QualifiedNameSyntaxError("A qualified name's option key part must not be empty");
            }

            // KEY := ['a'-'z', 'A'-'Z', '0'-'9', '-']
            static const char * valid_option_key_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789"
                    "-";

            auto pos = key.find_first_not_of(valid_option_key_characters);

            if (std::string::npos != pos)
            {
                throw QualifiedNameSyntaxError("'" + key + "' is not a valid option key part: Character '" + key[pos] + "' may not be used");
            }
        }

        OptionValue::OptionValue(const std::string & value) :
            _value(value)
        {
            if (value.empty())
            {
                throw QualifiedNameSyntaxError("A qualified name's option value part must not be empty");
            }

            // VALUE := ['a'-'z', 'A'-'Z', '0'-'9', '+', '-', '/', '.', '^', '_']
            static const char * valid_option_value_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789"
                    "+-/.^_";

            auto pos = value.find_first_not_of(valid_option_value_characters);

            if (std::string::npos != pos)
            {
                throw QualifiedNameSyntaxError("'" + value + "' is not a valid option value part: Character '" + value[pos] + "' may not be used");
            }
        }
    }

    QualifiedName::QualifiedName(const std::string & input) :
        _full(input),
        _prefix("null"),
        _name("empty"),
        _suffix(""),
        _options()
    {
        if (input.empty())
        {
            throw QualifiedNameSyntaxError("A qualified name must not be empty");
        }

        const auto pos_scope     = input.find("::");
        const auto pos_at        = input.find('@');
        const auto pos_semicolon = input.find(';');

        if (std::string::npos == pos_scope)
        {
            throw QualifiedNameSyntaxError("'" + input + "' is not a valid qualified name: No scope delimiter ('::') found");
        }

        // A valid prefix does not contain either a ';' or an '@'.
        _prefix = qnp::Prefix(input.substr(0, pos_scope));


        // Check that the suffix comes before the options list, prohibiting
        // e.g.:
        //
        //    foo::bar;key=value@suffix
        if ((std::string::npos != pos_at) && (std::string::npos != pos_semicolon) && (pos_at > pos_semicolon))
        {
            throw QualifiedNameSyntaxError("'" + input + "' is not a valid qualified name: The suffix part must preceed the options list part");
        }

        auto len_name = std::min(pos_at, pos_semicolon);
        if (std::string::npos != len_name)
        {
            len_name -= pos_scope + 2;
        }

        _name = qnp::Name(input.substr(pos_scope + 2, len_name));

        _str = _prefix.str() + "::" + _name.str();

        // The suffix is optional
        if (std::string::npos != pos_at)
        {
            auto len_suffix = pos_semicolon;
            if (std::string::npos != len_suffix)
            {
                len_suffix -= pos_at + 1;
            }

            _suffix = qnp::Suffix(input.substr(pos_at + 1, len_suffix));

            _str += "@" + _suffix.str();
        }

        auto pos_option_start = pos_semicolon;
        while (std::string::npos != pos_option_start)
        {
            auto pos_equal = input.find('=', pos_option_start + 1);
            if (std::string::npos == pos_equal)
            {
                throw QualifiedNameSyntaxError("'" + input + "' is not a valid qualified name: An option specification must contain a '='");
            }

            auto pos_next_comma = input.find(',', pos_equal + 1);

            qnp::OptionKey key(input.substr(pos_option_start + 1, pos_equal - pos_option_start - 1));
            qnp::OptionValue value(input.substr(pos_equal + 1, pos_next_comma - pos_equal - 1));

            _options.declare(key.str(), value.str());

            pos_option_start = pos_next_comma;
        }
    }

    QualifiedName::QualifiedName(const char * input) :
        QualifiedName(std::string(input))
    {
    }

    QualifiedName::QualifiedName(const QualifiedName & other) :
        _str(other._str),
        _full(other._full),
        _prefix(other._prefix),
        _name(other._name),
        _suffix(other._suffix),
        _options(other._options)
    {
    }

    QualifiedName::QualifiedName(const qnp::Prefix & p, const qnp::Name & n, const qnp::Suffix & s) :
        _str(p.str() + "::" + n.str() + (s.empty() ? std::string() : "@" + s.str())),
        _full(_str),
        _prefix(p),
        _name(n),
        _suffix(s),
        _options()
    {
    }

    QualifiedName::~QualifiedName()
    {
    }

    QualifiedNameSyntaxError::QualifiedNameSyntaxError(const std::string & msg) :
        Exception(msg)
    {
    }
}
