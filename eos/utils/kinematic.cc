/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2026 Danny van Dyk
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

#include <eos/utils/kinematic.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <limits>
#include <map>
#include <set>
#include <vector>

namespace eos
{
    template <> struct WrappedForwardIteratorTraits<Kinematics::KinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::vector<KinematicVariable>::iterator;
    };
    template class WrappedForwardIterator<Kinematics::KinematicVariableIteratorTag, const KinematicVariable>;

    template <> struct Implementation<Kinematics>
    {
            std::vector<double> variables_data;

            // Non-null only for derived variables.
            std::vector<std::shared_ptr<const KinematicHook>> hooks;

            std::map<std::string, unsigned> variables_map;

            std::vector<std::string> variables_names;

            std::map<std::string, unsigned> alias_map;

            std::vector<KinematicVariable> variables;

            double
            value(const unsigned & index) const
            {
                if (hooks[index]) [[unlikely]]
                {
                    return hooks[index]->evaluate();
                }

                return variables_data[index];
            }

            void
            set(const unsigned & index, const double & value)
            {
                if (hooks[index]) [[unlikely]]
                {
                    throw InternalError("Kinematic variable '" + variables_names[index] + "' is derived and cannot be set");
                }

                variables_data[index] = value;
            }
    };

    Kinematics::Kinematics() :
        PrivateImplementationPattern<Kinematics>(new Implementation<Kinematics>)
    {
    }

    Kinematics::Kinematics(const std::initializer_list<std::pair<std::string, double>> & variables) :
        PrivateImplementationPattern<Kinematics>(new Implementation<Kinematics>)
    {
        for (const auto & v : variables)
        {
            auto i(_imp->variables_map.find(v.first));

            if (_imp->variables_map.end() == i)
            {
                int index                    = _imp->variables_data.size();
                _imp->variables_map[v.first] = index;
                _imp->variables_data.push_back(v.second);
                _imp->hooks.push_back(nullptr);
                _imp->variables_names.push_back(v.first);
                _imp->variables.push_back(KinematicVariable(_imp, index, false));
            }
            else
            {
                _imp->variables_data[i->second]  = v.second;
                _imp->variables_names[i->second] = v.first;
            }
        }
    }

    Kinematics::~Kinematics() {}

    Kinematics
    Kinematics::clone() const
    {
        Kinematics result;

        result._imp.reset(new Implementation<Kinematics>(*_imp));

        for (auto & hook : result._imp->hooks)
        {
            if (hook)
            {
                hook = hook->clone(result);
            }
        }

        return result;
    }

    Kinematics
    Kinematics::operator+ (const Kinematics & rhs) const
    {
        Kinematics result = this->clone();

        for (unsigned i = 0; i < rhs._imp->variables_names.size(); ++i)
        {
            if (! rhs._imp->hooks[i])
            {
                result.declare(rhs._imp->variables_names[i], rhs._imp->variables_data[i]);
            }
        }

        // derived variables read the variables declared above
        for (unsigned i = 0; i < rhs._imp->variables_names.size(); ++i)
        {
            if (rhs._imp->hooks[i])
            {
                result.declare_derived(rhs._imp->variables_names[i], rhs._imp->hooks[i]->clone(result));
            }
        }

        for (const auto & x : rhs._imp->alias_map)
        {
            result.alias(x.first, rhs._imp->variables_names[x.second]);
        }

        return result;
    }

    bool
    Kinematics::operator== (const Kinematics & rhs) const
    {
        if (_imp->variables_map.size() != rhs._imp->variables_map.size())
        {
            return false;
        }

        for (auto l = _imp->variables_map.cbegin(), l_end = _imp->variables_map.cend(), r = rhs._imp->variables_map.cbegin(); l != l_end; ++l, ++r)
        {
            if (l->first != r->first)
            {
                return false;
            }

            if (_imp->value(l->second) != rhs._imp->value(r->second))
            {
                return false;
            }
        }

        if (_imp->alias_map.size() != rhs._imp->alias_map.size())
        {
            return false;
        }

        for (auto l = _imp->alias_map.cbegin(), l_end = _imp->alias_map.cend(), r = rhs._imp->alias_map.cbegin(); l != l_end; ++l, ++r)
        {
            if (l->first != r->first)
            {
                return false;
            }
            if (_imp->alias_map[l->first] != rhs._imp->alias_map[r->first])
            {
                return false;
            }
        }

        return true;
    }

    bool
    Kinematics::operator!= (const Kinematics & rhs) const
    {
        return ! (*this == rhs);
    }

    KinematicVariable
    Kinematics::operator[] (const std::string & name) const
    {
        // aliases take precedence over variables of the same name
        auto j(_imp->alias_map.find(name));

        if (_imp->alias_map.end() != j)
        {
            return KinematicVariable(_imp, j->second, true);
        }

        auto i(_imp->variables_map.find(name));

        if (_imp->variables_map.end() != i)
        {
            return KinematicVariable(_imp, i->second, false);
        }

        throw UnknownKinematicVariableError(name);
    }

    KinematicVariable
    Kinematics::variable(const std::string & name) const
    {
        auto i(_imp->variables_map.find(name));

        if (_imp->variables_map.end() == i)
        {
            throw UnknownKinematicVariableError(name);
        }

        return KinematicVariable(_imp, i->second, false);
    }

    bool
    Kinematics::has(const std::string & name) const
    {
        // aliases resolve through operator[] just as variables do, and must be found here too
        return (_imp->variables_map.cend() != _imp->variables_map.find(name)) || (_imp->alias_map.cend() != _imp->alias_map.find(name));
    }

    std::string
    Kinematics::as_string() const
    {
        std::string result;

        auto i(_imp->variables_map.cbegin()), i_end(_imp->variables_map.cend());
        if (i != i_end)
        {
            result = i->first + '=' + stringify(_imp->value(i->second));
            ++i;
        }

        for (; i != i_end; ++i)
        {
            result += ", " + i->first + '=' + stringify(_imp->value(i->second));
        }

        auto j(_imp->alias_map.cbegin()), j_end(_imp->alias_map.cend());
        for (; j != j_end; ++j)
        {
            result += ", " + j->first + "->" + _imp->variables_names[j->second];
        }


        return result;
    }

    void
    Kinematics::alias(const std::string & alias, const std::string & name)
    {
        const auto i(_imp->variables_map.find(name));
        const auto j(_imp->alias_map.find(name));
        bool       name_exists  = (_imp->variables_map.end() != i);
        bool       alias_exists = (_imp->alias_map.end() != j);

        if (! name_exists)
        {
            throw UnknownKinematicVariableError(name);
        }

        if (alias_exists)
        {
            throw DuplicateKinematicAliasError(alias, name);
        }

        _imp->alias_map[alias] = i->second;
    }

    void
    Kinematics::remove_alias(const std::string & alias)
    {
        const auto i(_imp->alias_map.find(alias));

        if (_imp->alias_map.end() == i)
        {
            throw UnknownKinematicAliasError(alias);
        }

        _imp->alias_map.erase(i);
    }

    void
    Kinematics::clear_aliases()
    {
        _imp->alias_map.clear();
    }

    KinematicVariable
    Kinematics::declare(const std::string & name, const double & value)
    {
        const auto i(_imp->variables_map.find(name));
        const auto j(_imp->alias_map.find(name));
        bool       regular    = (_imp->variables_map.end() != i);
        bool       alias      = (_imp->alias_map.end() != j);
        bool       undeclared = (! regular) && (! alias);

        if (undeclared)
        {
            int index                 = _imp->variables_data.size();
            _imp->variables_map[name] = index;
            _imp->variables_data.push_back(value);
            _imp->hooks.push_back(nullptr);
            _imp->variables_names.push_back(name);
            _imp->variables.push_back(KinematicVariable(_imp, index, false));

            return KinematicVariable(_imp, index, false);
        }
        else if (alias)
        {
            _imp->set(j->second, value);

            return KinematicVariable(_imp, j->second, true);
        }
        else // regular
        {
            _imp->set(i->second, value);

            return KinematicVariable(_imp, i->second, false);
        }
    }

    KinematicVariable
    Kinematics::declare_derived(const std::string & name, const std::shared_ptr<const KinematicHook> & hook)
    {
        if (this->has(name))
        {
            throw InternalError("Cannot declare the derived kinematic variable '" + name + "', since a variable of that name already exists");
        }

        int index                 = _imp->variables_data.size();
        _imp->variables_map[name] = index;
        _imp->variables_data.push_back(std::numeric_limits<double>::quiet_NaN());
        _imp->hooks.push_back(hook);
        _imp->variables_names.push_back(name);
        _imp->variables.push_back(KinematicVariable(_imp, index, false));

        return KinematicVariable(_imp, index, false);
    }

    std::shared_ptr<const KinematicHook>
    Kinematics::hook(const std::string & name) const
    {
        return _imp->hooks[(*this)[name].id()];
    }

    void
    Kinematics::set(const std::string & name, const double & value)
    {
        const auto i(_imp->variables_map.find(name));
        const auto j(_imp->alias_map.find(name));
        bool       regular    = (_imp->variables_map.end() != i);
        bool       alias      = (_imp->alias_map.end() != j);
        bool       undeclared = (! regular) && (! alias);

        if (undeclared)
        {
            throw UnknownKinematicVariableError(name);
        }

        if (alias)
        {
            _imp->set(j->second, value);
        }
        else // regular
        {
            _imp->set(i->second, value);
        }
    }

    Kinematics::KinematicVariableIterator
    Kinematics::begin() const
    {
        return _imp->variables.begin();
    }

    Kinematics::KinematicVariableIterator
    Kinematics::end() const
    {
        return _imp->variables.end();
    }

    KinematicVariable::KinematicVariable(const std::shared_ptr<Implementation<Kinematics>> & imp, unsigned index, bool is_alias) :
        _imp(imp),
        _index(index),
        _is_alias(is_alias)
    {
    }

    KinematicVariable::~KinematicVariable() {}

    MutablePtr
    KinematicVariable::clone() const
    {
        return MutablePtr(new KinematicVariable(_imp, _index, _is_alias));
    }

    KinematicVariable::operator double () const
    {
        return _imp->value(_index);
    }

    double
    KinematicVariable::operator() () const
    {
        return _imp->value(_index);
    }

    double
    KinematicVariable::evaluate() const
    {
        return _imp->value(_index);
    }

    const KinematicVariable &
    KinematicVariable::operator= (const double & value)
    {
        _imp->set(_index, value);

        return *this;
    }

    void
    KinematicVariable::set(const double & value)
    {
        _imp->set(_index, value);
    }

    const std::string &
    KinematicVariable::name() const
    {
        return _imp->variables_names[_index];
    }

    KinematicVariable::Id
    KinematicVariable::id() const
    {
        return _index;
    }

    /* KinematicUser */

    template <> struct WrappedForwardIteratorTraits<KinematicUser::ConstIteratorTag>
    {
            using UnderlyingIterator = std::set<KinematicVariable::Id>::const_iterator;
    };
    template class WrappedForwardIterator<KinematicUser::ConstIteratorTag, const KinematicVariable::Id>;

    KinematicUser::ConstIterator
    KinematicUser::begin_kinematics() const
    {
        return ConstIterator(_ids.cbegin());
    }

    KinematicUser::ConstIterator
    KinematicUser::end_kinematics() const
    {
        return ConstIterator(_ids.cend());
    }

    void
    KinematicUser::drop(const KinematicVariable::Id & id)
    {
        _ids.erase(id);
    }

    void
    KinematicUser::uses_kinematic(const KinematicVariable::Id & id)
    {
        _ids.insert(id);
    }

    void
    KinematicUser::uses_kinematic(const KinematicUser & other)
    {
        _ids.insert(other._ids.cbegin(), other._ids.cend());
    }

    UsedKinematicVariable::UsedKinematicVariable(const KinematicVariable & variable, KinematicUser & user) :
        KinematicVariable(variable)
    {
        user.uses_kinematic(variable.id());
    }

    UnknownKinematicVariableError::UnknownKinematicVariableError(const std::string & variable) throw() :
        Exception("Unknown kinematic variable: '" + variable + "'")
    {
    }

    DuplicateKinematicAliasError::DuplicateKinematicAliasError(const std::string & alias, const std::string & variable) throw() :
        Exception("Alias: '" + alias + "' cannot be used for variable: '" + variable + "' since it was already defined")
    {
    }

    UnknownKinematicAliasError::UnknownKinematicAliasError(const std::string & alias) throw() :
        Exception("Unknown kinematic alias: '" + alias + "'")
    {
    }
} // namespace eos
