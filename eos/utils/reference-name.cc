/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/utils/reference-name.hh>

namespace eos
{
    namespace rnp
    {
        Name::Name(const std::string & name) :
            _name(name)
        {
            if (name.empty())
            {
                throw ReferenceNameSyntaxError("A reference name's name part must not be empty");
            }

            static const std::string valid_characters =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789";

            static const std::string valid_characters_and_digits =
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ"
                    "0123456789";

            auto pos = valid_characters.find_first_of(name[0]);
            if (std::string::npos == pos)
            {
                throw ReferenceNameSyntaxError("'" + name + "' is not a valid name part: first character '" + name[0] + "' may not be used");
            }

            pos = name.find_first_not_of(valid_characters_and_digits);
            if (std::string::npos != pos)
            {
                throw ReferenceNameSyntaxError("'" + name + "' is not a valid name part: character '" + name[pos] + "' may not be used");
            }
        }

        Year::Year(const std::string & year) :
            _year(year)
        {
            if (year.empty())
            {
                throw ReferenceNameSyntaxError("A qualified name's year part must not be empty");
            }

            if (4 != year.length())
            {
                throw ReferenceNameSyntaxError("A qualified name's year part must be exactly 4 digits long");
            }

            // YEAR          := ['0'-'9'] ['0'-'9'] ['0'-'9'] ['0'-'9']
            static const char * valid_digits =
                    "0123456789";

            auto pos = year.find_first_not_of(valid_digits);
            if (std::string::npos != pos)
            {
                throw ReferenceNameSyntaxError("'" + year + "' is not a valid year part: character '" + year[pos] + "' may not be used");
            }
        }

        Index::Index(const std::string & index) :
            _index(index)
        {
            static const char * valid_index_characters =
                    "ABCDEFGHIJKLMNOPQRTSUVWXYZ";

            auto pos = index.find_first_not_of(valid_index_characters);

            if (index.empty())
            {
                throw ReferenceNameSyntaxError("A qualified name's index part must not be empty");
            }

            if (std::string::npos != pos)
            {
                throw ReferenceNameSyntaxError("'" + index + "' is not a valid index part: Character '" + index[pos] + "' may not be used");
            }
        }
    }

    ReferenceName::ReferenceName(const std::string & input) :
        _str(input),
        _name("A"),
        _year("0000"),
        _index("A")
    {
        if (input.empty())
        {
            throw ReferenceNameSyntaxError("A reference name must not be empty");
        }

        if (input.length() < 8)
        {
            throw ReferenceNameSyntaxError("A reference name must at least 8 characters long");
        }

        const auto pos_colon = input.find(':');
        if (std::string::npos == pos_colon)
        {
            throw ReferenceNameSyntaxError("A reference name must contain at least one ':'");
        }
        const auto pos_extra_colon = input.find(':', pos_colon + 1);
        if (std::string::npos != pos_extra_colon)
        {
            throw ReferenceNameSyntaxError("A reference name must contain exactly one ':'");
        }

        const auto pos_dash = input.find('-');
        if (std::string::npos == pos_dash)
        {
            throw ReferenceNameSyntaxError("A reference name must contain at least one '-'");
        }
        const auto pos_extra_dash = input.find('-', pos_dash + 1);
        if (std::string::npos != pos_extra_dash)
        {
            throw ReferenceNameSyntaxError("A reference name must contain exactly one '-'");
        }

        if (pos_dash < pos_colon)
        {
            throw ReferenceNameSyntaxError("In a reference name, the ':' must preceed the '-'");
        }

        _name  = rnp::Name(input.substr(0, pos_colon));
        _year  = rnp::Year(input.substr(pos_colon + 1, pos_dash - pos_colon - 1));
        _index = rnp::Index(input.substr(pos_dash + 1));
    }

    ReferenceName::ReferenceName(const char * input) :
        ReferenceName(std::string(input))
    {
    }

    ReferenceName::ReferenceName(const ReferenceName & other) :
        _str(other._str),
        _name(other._name),
        _year(other._year),
        _index(other._index)
    {
    }

    ReferenceName::~ReferenceName()
    {
    }

    ReferenceNameSyntaxError::ReferenceNameSyntaxError(const std::string & msg) :
        Exception(msg)
    {
    }
}
