/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2025      Florian Herren
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
#include <eos/utils/qualified-name-parts.hh>
#include <eos/utils/quantum-numbers.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <string>
#include <variant>
#include <vector>

namespace eos
{
    /*!
     * UnknownOptionError is thrown when an Options object does not contain a value for a given option key.
     */
    struct UnknownOptionError : public Exception
    {
            UnknownOptionError(const qnp::OptionKey & key) throw();
    };

    /*!
     * InvalidOptionValueError is thrown when the value passed to a known option is invalid.
     */
    struct InvalidOptionValueError : public Exception
    {
            InvalidOptionValueError(const qnp::OptionKey & key, const std::string & value, const std::string & allowed = "") throw();
    };

    /*!
     * UnspecifiedOptionError is thrown by an observable provider or similar when a mandatory option is not specified.
     */
    struct UnspecifiedOptionError : public Exception
    {
            UnspecifiedOptionError(const qnp::OptionKey & key, const std::string & allowed = "") throw();
    };

    /*!
     * Options keeps the set of all string options for any Observable.
     */
    class Options : public PrivateImplementationPattern<Options>
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
            Options(const std::initializer_list<std::pair<qnp::OptionKey, std::string>> & options);

            /// Destructor.
            ~Options();

            /// Equality comparison operator.
            bool operator== (const Options & rhs) const;

            /// Inequality comparison operator.
            bool operator!= (const Options & rhs) const;
            ///@}

            ///@name Access
            ///@{
            const std::string & operator[] (const qnp::OptionKey & key) const;

            bool has(const qnp::OptionKey & key) const;

            void declare(const qnp::OptionKey & key, const std::string & value = "");

            std::string get(const qnp::OptionKey & key, const std::string & default_value = "") const;

            std::string as_string() const;

            bool empty() const;
            ///@}

            ///@name Iteration over our options
            ///@{
            struct OptionIteratorTag;
            using OptionIterator = WrappedForwardIterator<OptionIteratorTag, const std::pair<const qnp::OptionKey, std::string>>;

            OptionIterator begin() const;
            OptionIterator end() const;
            ///@}
    };

    extern template class WrappedForwardIterator<Options::OptionIteratorTag, const std::pair<const qnp::OptionKey, std::string>>;

    /// Merge operator.
    Options operator+ (const Options & lhs, const Options & rhs);

    /*!
     * Metadata of an Option, providing key, allowed values, and default value.
     */
    struct OptionSpecification
    {
            OptionSpecification(const OptionSpecification &);
            OptionSpecification(const qnp::OptionKey & key_in, const std::vector<std::string> & allowed_values_in);
            OptionSpecification(const qnp::OptionKey & key_in, const std::vector<std::string> & allowed_values_in, const std::string & default_value_in);
            OptionSpecification(const qnp::OptionKey & key_in, const std::string & allowed_value_in);
            OptionSpecification(const qnp::OptionKey & key_in, const std::string & allowed_value_in, const std::string & default_value_in);

            ~OptionSpecification();

            OptionSpecification & operator= (const OptionSpecification &);

            using AllowedValues = std::variant<std::string, std::vector<std::string>>;

            qnp::OptionKey key;
            AllowedValues  allowed_values;
            std::string    default_value;
    };

    class SpecifiedOption
    {
        private:
            SpecifiedOption & operator= (const SpecifiedOption &);

        protected:
            OptionSpecification _specification;
            std::string         _value;

        public:
            SpecifiedOption(const SpecifiedOption &);
            SpecifiedOption(const Options & options, const OptionSpecification & specification);
            SpecifiedOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key);
            ~SpecifiedOption();

            const std::string & value() const;
    };

    class RestrictedOption : public SpecifiedOption
    {
        public:
            RestrictedOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key);
            ~RestrictedOption();
    };

    class BooleanOption : public SpecifiedOption
    {
        private:
            bool boolean_value;

        public:
            BooleanOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key = "true"_ok);
            ~BooleanOption();

            bool                value() const;
            const std::string & str() const;
    };

    class IntegerOption : public SpecifiedOption
    {
        private:
            double _int_value;

        public:
            IntegerOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key);
            ~IntegerOption();

            double              value() const;
            const std::string & str() const;
    };

    class FloatOption : public SpecifiedOption
    {
        private:
            double _float_value;

        public:
            FloatOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key);
            ~FloatOption();

            double              value() const;
            const std::string & str() const;
    };

    class LeptonFlavorOption : public RestrictedOption
    {
        public:
            LeptonFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key = "l"_ok);
            ~LeptonFlavorOption();

            LeptonFlavor        value() const;
            const std::string & str() const;
    };

    class QuarkFlavorOption : public RestrictedOption
    {
        public:
            QuarkFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key = "q"_ok);
            ~QuarkFlavorOption();

            QuarkFlavor         value() const;
            const std::string & str() const;
    };

    class LightMesonOption : public RestrictedOption
    {
        public:
            LightMesonOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key);
            ~LightMesonOption();

            LightMeson          value() const;
            const std::string & str() const;
    };

    class IsospinOption : public SpecifiedOption
    {
        private:
            Isospin _isospin_value;

        public:
            IsospinOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key = "I"_ok);
            ~IsospinOption();

            Isospin             value() const;
            const std::string & str() const;
    };

    class PartialWaveOption : public SpecifiedOption
    {
        private:
            PartialWave _partial_wave_value;

        public:
            PartialWaveOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key = "L"_ok);
            ~PartialWaveOption();

            PartialWave         value() const;
            const std::string & str() const;
    };
} // namespace eos

#endif
