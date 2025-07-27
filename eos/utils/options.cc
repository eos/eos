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

#include <eos/utils/destringify.hh>
#include <eos/utils/join.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <map>

namespace eos
{
    template <> struct WrappedForwardIteratorTraits<Options::OptionIteratorTag>
    {
            using UnderlyingIterator = std::map<qnp::OptionKey, std::string>::const_iterator;
    };
    template class WrappedForwardIterator<Options::OptionIteratorTag, const std::pair<const qnp::OptionKey, std::string>>;

    template <> struct Implementation<Options>
    {
            std::map<qnp::OptionKey, std::string> options;

            Implementation() {}

            Implementation(const std::initializer_list<std::pair<qnp::OptionKey, std::string>> & _options)
            {
                for (const auto & _option : _options)
                {
                    options.insert(_option);
                }
            }
    };

    Options::Options() :
        PrivateImplementationPattern<Options>(new Implementation<Options>)
    {
    }

    Options::Options(const std::initializer_list<std::pair<qnp::OptionKey, std::string>> & options) :
        PrivateImplementationPattern<Options>(new Implementation<Options>(options))
    {
    }

    Options::~Options() {}

    bool
    Options::operator== (const Options & rhs) const
    {
        if (_imp->options.size() != rhs._imp->options.size())
        {
            return false;
        }

        for (auto l = _imp->options.cbegin(), l_end = _imp->options.cend(), r = rhs._imp->options.cbegin(); l != l_end; ++l, ++r)
        {
            if (*l != *r)
            {
                return false;
            }
        }

        return true;
    }

    bool
    Options::operator!= (const Options & rhs) const
    {
        return ! (*this == rhs);
    }

    const std::string &
    Options::operator[] (const qnp::OptionKey & key) const
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() == i)
        {
            throw UnknownOptionError(key);
        }

        return i->second;
    }

    bool
    Options::has(const qnp::OptionKey & key) const
    {
        return _imp->options.end() != _imp->options.find(key);
    }

    void
    Options::declare(const qnp::OptionKey & key, const std::string & value)
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() != i)
        {
            i->second = value;
        }
        else
        {
            _imp->options[key] = value;
        }
    }

    std::string
    Options::get(const qnp::OptionKey & key, const std::string & default_value) const
    {
        auto i(_imp->options.find(key));

        if (_imp->options.end() == i)
        {
            return default_value;
        }

        return i->second;
    }

    std::string
    Options::as_string() const
    {
        std::string result;

        auto i(_imp->options.cbegin()), i_end(_imp->options.cend());

        if (i != i_end)
        {
            result += i->first.str() + '=' + i->second;
            ++i;
        }

        for (; i != i_end; ++i)
        {
            result += ',' + i->first.str() + '=' + i->second;
        }

        return result;
    }

    Options::OptionIterator
    Options::begin() const
    {
        return _imp->options.cbegin();
    }

    Options::OptionIterator
    Options::end() const
    {
        return _imp->options.cend();
    }

    bool
    Options::empty() const
    {
        return _imp->options.empty();
    }

    UnknownOptionError::UnknownOptionError(const qnp::OptionKey & key) throw() :
        Exception("Unknown option: '" + key.str() + "'")
    {
    }

    InvalidOptionValueError::InvalidOptionValueError(const qnp::OptionKey & key, const std::string & value, const std::string & allowed) throw() :
        Exception("Invalid value '" + value + "' for option: '" + key.str() + "'" + (allowed.empty() ? "" : ". Allowed values: '" + allowed + "'"))
    {
    }

    UnspecifiedOptionError::UnspecifiedOptionError(const qnp::OptionKey & key, const std::string & allowed) throw() :
        Exception("Mandatory option '" + key.str() + "' not specified'" + (allowed.empty() ? "" : ". Allowed values: '" + allowed + "'"))
    {
    }

    Options
    operator+ (const Options & lhs, const Options & rhs)
    {
        Options result;

        result._imp->options.insert(lhs._imp->options.cbegin(), lhs._imp->options.cend());

        /*
         * merge all options from rhs into result. Make sure to overwrite
         * existing lhs options with the same key.
         */
        for (auto o : rhs._imp->options)
        {
            auto retval = result._imp->options.insert(o);

            if (! retval.second)
            {
                retval.first->second = o.second;
            }
        }

        return result;
    }

    OptionSpecification::OptionSpecification(const OptionSpecification &) = default;

    OptionSpecification::OptionSpecification(const qnp::OptionKey & key_in, const std::vector<std::string> & allowed_values_in) :
        key(key_in),
        allowed_values(allowed_values_in)
    {
    }

    OptionSpecification::OptionSpecification(const qnp::OptionKey & key_in, const std::vector<std::string> & allowed_values_in, const std::string & default_value_in) :
        key(key_in),
        allowed_values(allowed_values_in),
        default_value(default_value_in)
    {
    }

    SpecifiedOption::SpecifiedOption(const SpecifiedOption &) = default;

    SpecifiedOption::SpecifiedOption(const Options & options, const OptionSpecification & specification) :
        _specification(specification)
    {
        if (! options.has(specification.key))
        {
            if ("" == specification.default_value)
            {
                throw UnspecifiedOptionError(specification.key, join(specification.allowed_values.cbegin(), specification.allowed_values.cend()));
            }

            _value = specification.default_value;
        }
        else
        {
            _value = options[specification.key];
        }
    }

    SpecifiedOption::SpecifiedOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        _specification("dummy"_ok, {}, "")
    {
        const auto s = std::find_if(specifications.cbegin(), specifications.cend(), [&](const auto & e) -> bool { return e.key == key; });

        if (specifications.cend() == s)
        {
            throw UnspecifiedOptionError("Options key '" + key.str() + "' is not specified in the options specifications");
        }

        *this = SpecifiedOption(options, *s);
    }

    SpecifiedOption::~SpecifiedOption() = default;

    SpecifiedOption & SpecifiedOption::operator= (const SpecifiedOption & other) = default;

    const std::string &
    SpecifiedOption::value() const
    {
        return _value;
    }

    RestrictedOption::RestrictedOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key)
    {
        if (std::find(_specification.allowed_values.cbegin(), _specification.allowed_values.cend(), _value) == _specification.allowed_values.cend())
        {
            throw InvalidOptionValueError(_specification.key, _value, join(_specification.allowed_values.cbegin(), _specification.allowed_values.cend()));
        }
    }

    RestrictedOption::~RestrictedOption() = default;

    BooleanOption::BooleanOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key),
        boolean_value(destringify<bool>(_value))
    {
    }

    BooleanOption::~BooleanOption() = default;

    bool
    BooleanOption::value() const
    {
        return boolean_value;
    }

    const std::string &
    BooleanOption::str() const
    {
        return _value;
    }

    IntegerOption::IntegerOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key),
        _int_value(destringify<int>(_value))
    {
    }

    IntegerOption::~IntegerOption() = default;

    double
    IntegerOption::value() const
    {
        return _int_value;
    }

    const std::string &
    IntegerOption::str() const
    {
        return _value;
    }

    FloatOption::FloatOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key),
        _float_value(destringify<double>(_value))
    {
    }

    FloatOption::~FloatOption() = default;

    double
    FloatOption::value() const
    {
        return _float_value;
    }

    const std::string &
    FloatOption::str() const
    {
        return _value;
    }

    LeptonFlavorOption::LeptonFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        RestrictedOption(options, specifications, key)
    {
    }

    LeptonFlavorOption::~LeptonFlavorOption() = default;

    LeptonFlavor
    LeptonFlavorOption::value() const
    {
        static const std::map<std::string, LeptonFlavor> map{
            {   "e", LeptonFlavor::electron },
            {  "mu",     LeptonFlavor::muon },
            { "tau",    LeptonFlavor::tauon }
        };

        const auto i = map.find(_value);
        if (map.cend() == i)
        {
            throw InternalError("Invalid lepton flavor '" + _value + "' encountered in LeptonFlavorOption::value()");
        }

        return i->second;
    }

    const std::string &
    LeptonFlavorOption::str() const
    {
        return _value;
    }

    QuarkFlavorOption::QuarkFlavorOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        RestrictedOption(options, specifications, key)
    {
    }

    QuarkFlavorOption::~QuarkFlavorOption() = default;

    QuarkFlavor
    QuarkFlavorOption::value() const
    {
        static const std::map<std::string, QuarkFlavor> map{
            { "u",      QuarkFlavor::up },
            { "d",    QuarkFlavor::down },
            { "s", QuarkFlavor::strange },
            { "c",   QuarkFlavor::charm },
            { "b",  QuarkFlavor::bottom },
            { "t",     QuarkFlavor::top }
        };

        const auto i = map.find(_value);
        if (map.cend() == i)
        {
            throw InternalError("Invalid quark flavor '" + _value + "' encountered in QuarkFlavorOption::value()");
        }

        return i->second;
    }

    const std::string &
    QuarkFlavorOption::str() const
    {
        return _value;
    }

    LightMesonOption::LightMesonOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        RestrictedOption(options, specifications, key)
    {
    }

    LightMesonOption::~LightMesonOption() = default;

    LightMeson
    LightMesonOption::value() const
    {
        static const std::map<std::string, LightMeson> map{
            {      "pi^0",     LightMeson::pi0 },
            {      "pi^+",  LightMeson::piplus },
            {      "pi^-", LightMeson::piminus },
            {       "K_d",      LightMeson::K0 },
            {    "Kbar_d",   LightMeson::K0bar },
            {       "K_S",      LightMeson::KS },
            {       "K_u",   LightMeson::Kplus },
            {    "Kbar_u",  LightMeson::Kminus },
            {       "eta",     LightMeson::eta },
            { "eta_prime",    LightMeson::etap }
        };

        const auto i = map.find(_value);
        if (map.cend() == i)
        {
            throw InternalError("Invalid light meson '" + _value + "' encountered in LightMesonOption::value()");
        }

        return i->second;
    }

    const std::string &
    LightMesonOption::str() const
    {
        return _value;
    }

    IsospinOption::IsospinOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key),
        _isospin_value(destringify<Isospin>(_value))
    {
        if (((_isospin_value ^ destringify<Isospin>(this->_specification.allowed_values.front())) & _isospin_value) != Isospin::none)
        {
            throw InvalidOptionValueError(_specification.key, _value, join(_specification.allowed_values.cbegin(), _specification.allowed_values.cend()));
        }
    }

    IsospinOption::~IsospinOption() = default;

    Isospin
    IsospinOption::value() const
    {
        return _isospin_value;
    }

    const std::string &
    IsospinOption::str() const
    {
        return _value;
    }

    PartialWaveOption::PartialWaveOption(const Options & options, const std::vector<OptionSpecification> & specifications, const qnp::OptionKey & key) :
        SpecifiedOption(options, specifications, key),
        _partial_wave_value(destringify<PartialWave>(_value))
    {
        if (((_partial_wave_value ^ destringify<PartialWave>(this->_specification.allowed_values.front())) & _partial_wave_value) != PartialWave::none)
        {
            throw InvalidOptionValueError(_specification.key, _value, join(_specification.allowed_values.cbegin(), _specification.allowed_values.cend()));
        }
    }

    PartialWaveOption::~PartialWaveOption() = default;

    PartialWave
    PartialWaveOption::value() const
    {
        return _partial_wave_value;
    }

    const std::string &
    PartialWaveOption::str() const
    {
        return _value;
    }
} // namespace eos
