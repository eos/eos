/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2024 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_OPTIONS_HH
#define EOS_GUARD_EOS_UTILS_OPTIONS_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/quantum-numbers.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <string>
#include <vector>

namespace eos
{
    /*!
     * UnknownOptionError is thrown when an Options object does not contain a value for a given option key.
     */
    struct UnknownOptionError :
        public Exception
    {
        UnknownOptionError(const std::string & key) throw ();
    };

    /*!
     * InvalidOptionValueError is thrown when the value passed to a known option is invalid.
     */
    struct InvalidOptionValueError :
        public Exception
    {
        InvalidOptionValueError(const std::string & key, const std::string & value, const std::string & allowed = "") throw ();
    };

    /*!
     * UnspecifiedOptionError is thrown by an observable provider or similar when a mandatory option is not specified.
     */
    struct UnspecifiedOptionError :
        public Exception
    {
        UnspecifiedOptionError(const std::string & key, const std::string & allowed = "") throw ();
    };

    /*!
     * Options keeps the set of all string options for any Observable.
     */
    class Options :
        public PrivateImplementationPattern<Options>
    {
        public:
            friend Options operator+ (const Options &, const Options &);

            ///@name Basic Functions
            ///@{

            /// Constructor.
            Options();

            /*!
             * Constructor.
             *
             * Create an instance of Kinematics with a given set of initial
             * options.
             *
             * @param options The set of initial options from which this object shall be constructed.
             */
            Options(const std::initializer_list<std::pair<std::string, std::string>> & options);

            /// Destructor.
            ~Options();

            /// Equality comparison operator.
            bool operator== (const Options & rhs) const;

            /// Inequality comparison operator.
            bool operator!= (const Options & rhs) const;
            ///@}

            ///@name Access
            ///@{
            const std::string & operator[] (const std::string & key) const;

            bool has(const std::string & key) const;

            void declare(const std::string & key, const std::string & value = "");

            std::string get(const std::string & key, const std::string & default_value = "") const;

            std::string as_string() const;

            bool empty() const;
            ///@}

            ///@name Iteration over our options
            ///@{
            struct OptionIteratorTag;
            using OptionIterator = WrappedForwardIterator<OptionIteratorTag, const std::pair<const std::string, std::string>>;

            OptionIterator begin() const;
            OptionIterator end() const;
            ///@}
    };

    extern template class WrappedForwardIterator<Options::OptionIteratorTag, const std::pair<const std::string, std::string>>;

    /// Merge operator.
    Options operator+ (const Options & lhs, const Options & rhs);

    /*!
     * Metadata of an Option, providing key, allowed values, and default value.
     */
    struct OptionSpecification
    {
        std::string key;
        std::vector<std::string> allowed_values;
        std::string default_value;
    };

    class SpecifiedOption
    {
        protected:
            std::string _value;

        public:
            SpecifiedOption(const Options & options, const OptionSpecification & specification);
            SpecifiedOption(const Options & options, const std::vector<OptionSpecification> & specifications, const std::string & key);
            ~SpecifiedOption();

            const std::string & value() const;
    };

    class BooleanOption :
        public SpecifiedOption
    {
        private:
            bool boolean_value;

        public:
            BooleanOption(const Options & options, const std::vector<OptionSpecification> & specifications, const std::string & key = "true");
            ~BooleanOption();

            bool value() const;
            const std::string & str() const;
    };

    class FloatOption :
        public SpecifiedOption
    {
        private:
            double _float_value;

        public:
            FloatOption(const Options & options, const std::vector<OptionSpecification> & specifications, const std::string & key);
            ~FloatOption();

            double value() const;
            const std::string & str() const;
    };

    class LeptonFlavorOption :
        public SpecifiedOption
    {
        public:
            LeptonFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const std::string & key = "l");
            ~LeptonFlavorOption();

            LeptonFlavor value() const;
            const std::string & str() const;
    };

    class QuarkFlavorOption :
        public SpecifiedOption
    {
        public:
            QuarkFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const std::string & key = "q");
            ~QuarkFlavorOption();

            QuarkFlavor value() const;
            const std::string & str() const;
    };
}

#endif
