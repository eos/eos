/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
 * Copyright (c) 2013 Frederik Beaujean
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

#include <eos/statistics/simple-parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <map>

namespace eos
{
    template <>
    struct WrappedForwardIteratorTraits<SimpleParameters::IteratorTag>
    {
        typedef std::vector<ParameterDescription>::const_iterator UnderlyingIterator;
    };

    template <>
    struct Implementation<SimpleParameters>
    {
        // forbid parameters with same name
        std::map<std::string, SimpleParameter::Index> parameters_map;
        std::shared_ptr<std::vector<double>> values;
        std::vector<ParameterDescription> defs;

        Implementation() :
            values(new std::vector<double>)
        {

        }

        SimpleParameter & declare(const std::string & name, const double & min,
                                const double & max, bool nuisance=false)
        {
            auto i(parameters_map.find(name));

            if (parameters_map.end() == i)
            {
                SimpleParameter::Index id = defs.size();
                parameters_map[name] = id;
                values->push_back(0);
                SimpleParameter * p = new SimpleParameter(name, id, values);
                defs.push_back(ParameterDescription{ MutablePtr(p), min, max, nuisance });
                return *p;
            }
            else
            {
                return *(static_cast<SimpleParameter *>(defs[i->second].parameter.get()));
            }
        }
    };

    SimpleParameter::SimpleParameter(const std::string & name, const Index & index,
                                     const std::shared_ptr<std::vector<double>> & parameters) :
        _name(name),
        _index(index),
        _parameters(parameters)
    {
    }

    SimpleParameter::~SimpleParameter()
    {
    }

    MutablePtr
    SimpleParameter::clone() const
    {
        return MutablePtr(new SimpleParameter(_name, _index, _parameters));
    }

    SimpleParameter::operator double () const
    {
        return (*_parameters)[_index];
    }

    double
    SimpleParameter::operator() () const
    {
        return (*_parameters)[_index];
    }

    double
    SimpleParameter::evaluate() const
    {
        return (*_parameters)[_index];
    }

    const Mutable &
    SimpleParameter::operator= (const double & value)
    {
        (*_parameters)[_index] = value;

        return *this;
    }

    void
    SimpleParameter::set(const double & value)
    {
        (*_parameters)[_index] = value;
    }

    const std::string &
    SimpleParameter::name() const
    {
        return _name;
    }

    SimpleParameters::SimpleParameters() :
        PrivateImplementationPattern<SimpleParameters>(new Implementation<SimpleParameters>())
    {
    }

    SimpleParameters::~SimpleParameters()
    {
    }

    SimpleParameters
    SimpleParameters::clone() const
    {
        SimpleParameters result;

        // copy parameters
        for (auto & d : _imp->defs)
        {
            result.declare(d.parameter->name(), d.min, d.max, d.nuisance);
        }

        // copy values
        std::copy(_imp->values->begin(), _imp->values->end(), result._imp->values->begin());

        return result;
    }

    SimpleParameters::Iterator
    SimpleParameters::begin() const
    {
        return _imp->defs.cbegin();
    }

    SimpleParameters::Iterator
    SimpleParameters::end() const
    {
        return _imp->defs.cend();
    }

    SimpleParameter &
    SimpleParameters::declare(const std::string & name, const double & min,
                              const double & max, bool nuisance)
    {
        return _imp->declare(name, min, max, nuisance);
    }

    SimpleParameter &
    SimpleParameters::operator[] (const std::string & name) const
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        return *(static_cast<SimpleParameter *>(_imp->defs[i->second].parameter.get()));
    }

    SimpleParameter &
    SimpleParameters::operator[] (const SimpleParameter::Index & id) const
    {
        if (id >= _imp->defs.size())
            throw InternalError("Parameters::operator[] (Parameter::Id): invalid id '" + stringify(id) + "'");

        return *(static_cast<SimpleParameter *>(_imp->defs[id].parameter.get()));
    }

    bool
    SimpleParameters::operator!= (const SimpleParameters & rhs) const
    {
        return _imp->values.get() != rhs._imp->values.get();
    }

    const std::vector<double> &
    SimpleParameters::values() const
    {
        return *_imp->values;
    }

    template class WrappedForwardIterator<SimpleParameters::IteratorTag, const ParameterDescription>;
}
