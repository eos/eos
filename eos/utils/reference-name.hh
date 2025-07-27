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

#ifndef EOS_GUARD_EOS_UTILS_REFERENCE_HH
#define EOS_GUARD_EOS_UTILS_REFERENCE_HH 1

#include <eos/utils/exception.hh>

#include <string>

namespace eos
{
    namespace rnp
    {
        class Name
        {
            private:
                std::string _name;

            public:
                Name(const std::string &);

                const std::string &
                str() const
                {
                    return _name;
                }

                inline bool
                operator< (const Name & rhs) const
                {
                    return this->_name < rhs._name;
                }
        };

        class Year
        {
            private:
                std::string _year;

            public:
                Year(const std::string &);

                const std::string &
                str() const
                {
                    return _year;
                }

                inline bool
                operator< (const Year & rhs) const
                {
                    return this->_year < rhs._year;
                }
        };

        class Index
        {
            private:
                std::string _index;

            public:
                Index(const std::string &);

                const std::string &
                str() const
                {
                    return _index;
                }

                bool
                empty() const
                {
                    return _index.empty();
                }

                inline bool
                operator< (const Index & rhs) const
                {
                    return this->_index < rhs._index;
                }
        };
    } // namespace rnp

    /*
     * Holds a syntactically-correct reference name.
     *
     * The syntax of a reference name follows:
     *
     *   NAME ':' YEAR INDEX
     *
     * with
     *
     *   NAME  := ['a'-'z', 'A'-'Z'] ['a'-'z', 'A'-'Z', '0'-'9']*
     *   YEAR  := ['0'-'9'] ['0'-'9'] ['0'-'9'] ['0'-'9']
     *   INDEX := ['A'-'Z']+
     */
    class ReferenceName
    {
            friend std::ostream & operator<< (std::ostream &, const ReferenceName &);

        private:
            std::string _str;
            rnp::Name   _name;
            rnp::Year   _year;
            rnp::Index  _index;

        public:
            ReferenceName(const std::string & name);
            ReferenceName(const char * name);
            ReferenceName(const ReferenceName & other);
            ~ReferenceName();

            inline const std::string &
            str() const
            {
                return _str;
            }

            inline const rnp::Name &
            name_part() const
            {
                return _name;
            }

            inline const rnp::Year &
            year_part() const
            {
                return _year;
            }

            inline const rnp::Index &
            index_part() const
            {
                return _index;
            }

            inline bool
            operator< (const ReferenceName & rhs) const
            {
                return this->_str < rhs._str;
            }

            inline bool
            operator== (const ReferenceName & rhs) const
            {
                return this->_str == rhs._str;
            }

            inline bool
            operator!= (const ReferenceName & rhs) const
            {
                return this->_str != rhs._str;
            }
    };

    inline ReferenceName
    operator""_rn (const char * c, size_t)
    {
        return ReferenceName(c);
    }

    class ReferenceNameSyntaxError : public Exception
    {
        public:
            ReferenceNameSyntaxError(const std::string & msg);
    };

    inline std::ostream &
    operator<< (std::ostream & lhs, const ReferenceName & rhs)
    {
        lhs << rhs._str;

        return lhs;
    }
} // namespace eos

#endif
