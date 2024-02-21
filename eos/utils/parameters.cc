/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2024 Danny van Dyk
 * Copyright (c) 2021 Philip Lüghausen
 * Copyright (c) 2010 Christian Wacker
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

#include <config.h>

#include <eos/utils/cartesian-product.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <random>
#include <vector>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
#include <yaml-cpp/yaml.h>

#include <iostream>

namespace fs = boost::filesystem;

namespace eos
{
    bool
    operator== (const ParameterDescription & lhs, const ParameterDescription & rhs)
    {
    	if (lhs.min != rhs.min)
    		return false;
    	if (lhs.max != rhs.max)
    		return false;
    	if (lhs.nuisance!= rhs.nuisance)
    		return false;
    	if (lhs.parameter->name() != rhs.parameter->name())
    		return false;

    	return true;
    }

    /* ParameterGroup */

    template <>
    struct Implementation<ParameterGroup>
    {
        std::string name;

        std::string description;

        std::vector<Parameter> entries;

        Implementation(const std::string & name, const std::string & description, std::vector<Parameter> && entries) :
            name(name),
            description(description),
            entries(entries)
        {
        }
    };

    template <>
    struct WrappedForwardIteratorTraits<ParameterGroup::ParameterIteratorTag>
    {
        using UnderlyingIterator = std::vector<Parameter>::const_iterator;
    };
    template class WrappedForwardIterator<ParameterGroup::ParameterIteratorTag, const Parameter>;

    ParameterGroup::ParameterGroup(Implementation<ParameterGroup> * imp) :
        PrivateImplementationPattern<ParameterGroup>(imp)
    {
    }

    ParameterGroup::~ParameterGroup() = default;

    ParameterGroup::ParameterIterator
    ParameterGroup::begin() const
    {
        return _imp->entries.begin();
    }

    ParameterGroup::ParameterIterator
    ParameterGroup::end() const
    {
        return _imp->entries.end();
    }

    const std::string &
    ParameterGroup::name() const
    {
        return _imp->name;
    }

    const std::string &
    ParameterGroup::description() const
    {
        return _imp->description;
    }

    /* ParameterSection */

    template <>
    struct Implementation<ParameterSection>
    {
        std::string name;

        std::string description;

        std::vector<ParameterGroup> groups;

        Implementation(const std::string & name, const std::string & description, const std::vector<ParameterGroup> & groups) :
            name(name),
            description(description),
            groups(groups)
        {
        }
    };

    template <>
    struct WrappedForwardIteratorTraits<ParameterSection::GroupIteratorTag>
    {
        using UnderlyingIterator = std::vector<ParameterGroup>::const_iterator;
    };
    template class WrappedForwardIterator<ParameterSection::GroupIteratorTag, const ParameterGroup &>;

    ParameterSection::ParameterSection(Implementation<ParameterSection> * imp) :
        PrivateImplementationPattern<ParameterSection>(imp)
    {
    }

    ParameterSection::~ParameterSection() = default;

    ParameterSection::GroupIterator
    ParameterSection::begin() const
    {
        return _imp->groups.begin();
    }

    ParameterSection::GroupIterator
    ParameterSection::end() const
    {
        return _imp->groups.end();
    }

    const std::string &
    ParameterSection::name() const
    {
        return _imp->name;
    }

    const std::string &
    ParameterSection::description() const
    {
        return _imp->description;
    }

    struct Parameter::Template
    {
        QualifiedName name;

        double min, central, max;

        std::string latex;

        Unit unit;
    };

    struct Parameter::Data :
        Parameter::Template
    {
        double value, generator_value;

        Parameter::Id id;

        Data(const Parameter::Template & t, const Parameter::Id & i) :
            Parameter::Template(t),
            value(t.central),
            generator_value(0.0),
            id(i)
        {
        }
    };

    struct Parameters::Data
    {
        std::vector<Parameter::Data> data;
    };

    template <>
    struct WrappedForwardIteratorTraits<Parameters::IteratorTag>
    {
        using UnderlyingIterator = std::vector<Parameter>::iterator;
    };
    template class WrappedForwardIterator<Parameters::IteratorTag, Parameter>;

    template <>
    struct WrappedForwardIteratorTraits<Parameters::SectionIteratorTag>
    {
        using UnderlyingIterator = std::vector<ParameterSection>::const_iterator;
    };
    template class WrappedForwardIterator<Parameters::SectionIteratorTag, const ParameterSection &>;

    class ParameterDefaults :
        public InstantiationPolicy<ParameterDefaults, Singleton>
    {
        private:
            std::shared_ptr<Parameters::Data> _data;

            std::map<QualifiedName, unsigned> _map;

            std::vector<ParameterSection> _sections;

            ParameterDefaults() :
                _data(new Parameters::Data)
            {
                _load_defaults();
            }

            ~ParameterDefaults() = default;

            void
            _load_defaults()
            {
                Context ctx("When loading parameter defaults");

                fs::path base;
                if (std::getenv("EOS_TESTS_PARAMETERS"))
                {
                    std::string envvar = std::string(std::getenv("EOS_TESTS_PARAMETERS"));
                    base = fs::system_complete(envvar);
                }
                else if (std::getenv("EOS_HOME"))
                {
                    std::string envvar = std::string(std::getenv("EOS_HOME"));
                    base = fs::system_complete(envvar) / "parameters";
                }
                else
                {
                    base = fs::system_complete(EOS_DATADIR "/eos/parameters/");
                }

                if (! fs::exists(base))
                {
                    throw InternalError("Could not find the parameter input files, '" + base.string() + "' does not exist");
                }

                if (! fs::is_directory(base))
                {
                    throw InternalError("Expect '" + base.string() + " to be a directory");
                }

                unsigned idx = _data->data.size();
                for (fs::directory_iterator f(base), f_end ; f != f_end ; ++f)
                {
                    auto file_path = f->path();

                    if (! fs::is_regular_file(status(file_path)))
                        continue;

                    if (".yaml" != file_path.extension().string())
                        continue;

                    std::string file = file_path.string();
                    Context ctx("When parsing parameter file '" + file + "'");

                    try
                    {
                        YAML::Node root_node = YAML::LoadFile(file);
                        std::vector<ParameterGroup> section_groups;

                        // parse the section metadata
                        auto section_title_node = root_node["title"];
                        if (! section_title_node)
                            throw ParameterInputFileNodeError(file, "/", "has no entry named 'title'");

                        if (YAML::NodeType::Scalar != section_title_node.Type())
                            throw ParameterInputFileNodeError(file, "title", "is not a scalar");
                        std::string section_title = section_title_node.as<std::string>();

                        Context ctx("When parsing metadata for section '" + section_title + "'");

                        auto section_desc_node = root_node["description"];
                        if (! section_desc_node)
                            throw ParameterInputFileNodeError(file, "/", "has no entry named 'description'");
                        if (YAML::NodeType::Scalar != section_desc_node.Type())
                            throw ParameterInputFileNodeError(file, "description", "is not a scalar");
                        std::string section_desc = section_desc_node.as<std::string>();

                        auto section_groups_node = root_node["groups"];
                        if (! section_groups_node)
                            throw ParameterInputFileNodeError(file, "/", "has no entry named 'group'");
                        if (YAML::NodeType::Sequence != section_groups_node.Type())
                            throw ParameterInputFileNodeError(file, "groups", "is not a sequence");

                        // parse the section's groups
                        for (auto && group_node : section_groups_node)
                        {
                            std::vector<Parameter> group_parameters;

                            auto group_title_node = group_node["title"];
                            if (! group_title_node)
                                throw ParameterInputFileNodeError(file, "", "has no entry named 'title'");
                            if (YAML::NodeType::Scalar != group_title_node.Type())
                                throw ParameterInputFileNodeError(file, "title", "is not a scalar");
                            std::string group_title = group_title_node.as<std::string>();

                            Context ctx("When parsing metadata for group '" + group_title + "'");

                            auto group_desc_node = group_node["description"];
                            if (! group_desc_node)
                                throw ParameterInputFileNodeError(file, group_title, "has no entry named 'description'");
                            if (YAML::NodeType::Scalar != group_desc_node.Type())
                                throw ParameterInputFileNodeError(file, "'" + group_title + "'.description", "is not a scalar");
                            std::string group_desc = group_desc_node.as<std::string>();

                            auto group_parameters_node = group_node["parameters"];
                            if (! group_parameters_node)
                                throw ParameterInputFileNodeError(file, group_title, "has no entry named 'parameters'");
                            if (YAML::NodeType::Map != group_parameters_node.Type())
                                throw ParameterInputFileNodeError(file, "'" + group_title + "'.parameters", "is not a map");

                            // parse the group's parameters
                            for (auto && p : group_parameters_node)
                            {
                                std::string name = p.first.Scalar();
                                Context ctx("When parsing metadata for parameter '" + name + "'");

                                double central, min, max;

                                std::string latex;

                                auto central_node = p.second["central"];
                                if (! central_node)
                                    throw ParameterInputFileNodeError(file, name, "has no entry named 'central'");
                                if (YAML::NodeType::Scalar != central_node.Type())
                                    throw ParameterInputFileNodeError(file, name + ".central", "is not a scalar");
                                central = central_node.as<double>();

                                auto min_node = p.second["min"];
                                if (! min_node)
                                    throw ParameterInputFileNodeError(file, name, "has no entry named 'min'");
                                if (YAML::NodeType::Scalar != min_node.Type())
                                    throw ParameterInputFileNodeError(file, name + ".min", "is not a scalar");
                                min = min_node.as<double>();

                                auto max_node = p.second["max"];
                                if (! max_node)
                                    throw ParameterInputFileNodeError(file, name, "has no entry named 'max'");
                                if (YAML::NodeType::Scalar != max_node.Type())
                                    throw ParameterInputFileNodeError(file, name + ".max", "is not a scalar");
                                max = max_node.as<double>();

                                auto latex_node = p.second["latex"];

                                Unit unit = Unit::Undefined();
                                auto unit_node = p.second["unit"];
                                if (unit_node)
                                {
                                    if (YAML::NodeType::Scalar != unit_node.Type())
                                        throw ParameterInputFileNodeError(file, name + ".unit", "is not a scalar");

                                    unit = Unit(unit_node.as<std::string>());
                                }

                                if (name.find("%") == std::string::npos) // The parameter is not templated
                                {
                                    if (_map.end() != _map.find(name))
                                    {
                                        throw ParameterInputDuplicateError(file, name);
                                    }

                                    if (latex_node)
                                    {
                                        if (YAML::NodeType::Scalar != latex_node.Type())
                                            throw ParameterInputFileNodeError(file, name + ".latex", "is not a scalar");

                                        latex = latex_node.as<std::string>();
                                    }

                                    _data->data.push_back(Parameter::Data(Parameter::Template { QualifiedName(name), min, central, max, latex, unit }, idx));
                                    _map[name] = idx;
                                    group_parameters.push_back(Parameter(_data, idx));

                                    ++idx;
                                }
                                else // The parameter is templated
                                {
                                    Context ctx("When parsing metadata for a templated parameter");

                                    auto matrix_node = p.second["matrix"];
                                    if (! matrix_node)
                                    {
                                        throw ParameterInputFileNodeError(file, name, "is templated but doesn't have substitutions");
                                    }
                                    else
                                    {
                                        std::map<std::string, std::string> latex_map;
                                        if (latex_node)
                                        {
                                            if (YAML::NodeType::Scalar == latex_node.Type())
                                            {
                                                latex = latex_node.as<std::string>();
                                            }
                                            else if (YAML::NodeType::Sequence == matrix_node.Type())
                                            {
                                                    if (YAML::NodeType::Map != latex_node.Type())
                                                        throw ParameterInputFileNodeError(file, name + ".latex", "is not a map");

                                                    if ((! latex_node["template"]) || (! latex_node["map"]))
                                                    {
                                                        throw ParameterInputFileNodeError(file, name + ".latex", "is incomplete, needs 'template' and 'map' subkeys");
                                                    }
                                                    latex = latex_node["template"].as<std::string>();
                                                    latex_map = latex_node["map"].as<std::map<std::string, std::string>>();
                                            }
                                            else
                                                throw ParameterInputFileNodeError(file, name + ".matrix", "is not a scalar nor a sequence");
                                        }

                                        CartesianProduct<std::vector<std::string>> cp;

                                        // Parse the parameter substitutions and add them to the cartesian product
                                        for (auto && substitution : matrix_node)
                                        {
                                            std::vector<std::string> instances;

                                            for (auto && instance : substitution)
                                            {
                                                instances.push_back(instance.as<std::string>());
                                            }

                                            cp.over(instances);
                                        }

                                        for (auto cp_it = cp.begin() ; cp.end() != cp_it ; ++cp_it)
                                        {
                                            boost::format templated_name(name);
                                            boost::format templated_latex(latex);

                                            for (auto && i : *cp_it)
                                            {
                                                templated_name  = templated_name  % i;

                                                std::string mapped_i;
                                                auto latex_map_it = latex_map.find(i);
                                                if (latex_map.end() != latex_map_it)
                                                {
                                                    mapped_i = latex_map_it->second;
                                                }
                                                else
                                                {
                                                    mapped_i = i;
                                                }
                                                templated_latex = templated_latex % mapped_i;
                                            }

                                            QualifiedName qn(templated_name.str());

                                            if (_map.end() != _map.find(qn))
                                            {
                                                throw ParameterInputDuplicateError(file, qn.str());
                                            }

                                            _data->data.push_back(Parameter::Data(Parameter::Template { qn, min, central, max, templated_latex.str(), unit }, idx));
                                            _map[templated_name.str()] = idx;
                                            group_parameters.push_back(Parameter(_data, idx));

                                            ++idx;
                                        }
                                    }
                                }
                            }

                            section_groups.push_back(ParameterGroup(new Implementation<ParameterGroup>(group_title, group_desc, std::move(group_parameters))));
                        }
                        _sections.push_back(ParameterSection(new Implementation<ParameterSection>(section_title, section_desc, std::move(section_groups))));
                    }
                    catch (std::exception & e)
                    {
                        throw ParameterInputFileParseError(file, e.what());
                    }
                }
            }

        public:
            friend class InstantiationPolicy<ParameterDefaults, Singleton>;

            const Parameters::Data & data() const
            {
                return *_data;
            }

            const std::map<QualifiedName, unsigned> & map() const
            {
                return _map;
            }

            const std::vector<ParameterSection> & sections() const
            {
                return _sections;
            }

            Parameter::Id declare(const QualifiedName & key, const Parameter::Template & value)
            {
                unsigned idx = _data->data.size();
                _data->data.push_back(Parameter::Data{ value, idx });
                _map[key] = idx;

                return idx;
            }

            void redirect(const QualifiedName & name, const Parameter::Id & id)
            {
                auto i = _map.find(name);
                if (_map.end() == i)
                    throw UnknownParameterError(name);

                i->second = id;
            }
    };

    template <>
    struct Implementation<Parameters>
    {
        std::shared_ptr<Parameters::Data> parameters_data;

        std::map<QualifiedName, unsigned> parameters_map;

        std::vector<Parameter> parameters;

        Implementation() :
            parameters_data(new Parameters::Data(ParameterDefaults::instance()->data())),
            parameters_map(ParameterDefaults::instance()->map())
        {
            for (auto i = parameters_data->data.begin(), i_end = parameters_data->data.end() ; i != i_end ; ++i)
            {
                parameters.push_back(Parameter(parameters_data, i->id));
            }
        }

        Implementation(const Implementation & other) :
            parameters_data(new Parameters::Data(*other.parameters_data)),
            parameters_map(other.parameters_map)
        {
            parameters.reserve(other.parameters.size());
            for (unsigned i = 0 ; i != parameters.size() ; ++i)
            {
                parameters.push_back(Parameter(parameters_data, i));
            }
        }

        void
        override_from_file(const std::string & file)
        {
            fs::path file_path(file);

            if (! fs::is_regular_file(fs::status(file_path)))
            {
                throw ParameterInputFileParseError(file, "expect the parameter file to be a regular file");
            }

            try
            {
                YAML::Node node = YAML::LoadFile(file);

                for (auto && p : node)
                {
                    std::string name = p.first.Scalar();

                    if ("@metadata@" == name)
                        continue;

                    double central, min = -std::numeric_limits<double>::max(), max = +std::numeric_limits<double>::max();

                    std::string latex;

                    Unit unit = Unit::Undefined();

                    bool has_min = false, has_max = false, has_latex = false, has_unit = false;

                    if (! p.second["central"])
                    {
                        throw ParameterInputFileNodeError(file, name, "has no entry named 'central'");
                    }
                    else if (YAML::NodeType::Scalar != p.second["central"].Type())
                    {
                        throw ParameterInputFileNodeError(file, name + ".central", "is not a scalar");
                    }
                    central = p.second["central"].as<double>();

                    if (p.second["min"])
                    {
                        if (YAML::NodeType::Scalar != p.second["min"].Type())
                        {
                            throw ParameterInputFileNodeError(file, name + ".min", "is not a scalar");
                        }
                        min = p.second["min"].as<double>();
                        has_min = true;
                    }

                    if (p.second["max"])
                    {
                        if (YAML::NodeType::Scalar != p.second["max"].Type())
                        {
                            throw ParameterInputFileNodeError(file, name + ".max", "is not a scalar");
                        }
                        max = p.second["max"].as<double>();
                        has_max = true;
                    }

                    if (p.second["latex"])
                    {
                        if (YAML::NodeType::Scalar != p.second["latex"].Type())
                            throw ParameterInputFileNodeError(file, name + ".latex", "is not a scalar");

                        latex = p.second["latex"].as<std::string>();
                        has_latex = true;
                    }

                    auto unit_node = p.second["unit"];
                    if (unit_node)
                    {
                        if (YAML::NodeType::Scalar != unit_node.Type())
                            throw ParameterInputFileNodeError(file, name + ".unit", "is not a scalar");

                        unit = Unit(unit_node.as<std::string>());
                    }

                    auto i = parameters_map.find(name);
                    if (parameters_map.end() != i)
                    {
                        Log::instance()->message("[parameters.override]", ll_informational)
                            << "Overriding existing parameter '" << name << "' with central value '" << central << "'";

                        parameters_data->data[i->second].value = central;
                        if (has_min)
                        {
                            parameters_data->data[i->second].min = min;
                        }
                        if (has_max)
                        {
                            parameters_data->data[i->second].max = max;
                        }
                        if (has_latex)
                        {
                            parameters_data->data[i->second].latex = latex;
                        }
                        if (has_unit)
                        {
                            parameters_data->data[i->second].unit = unit;
                        }
                    }
                    else
                    {
                        Log::instance()->message("[parameters.override]", ll_informational)
                            << "Adding new parameter '" << name << "' with central value '" << central << "'";

                        if (! has_min)
                        {
                            min = central;
                        }
                        if (! has_max)
                        {
                            max = central;
                        }

                        auto idx = parameters_data->data.size();
                        parameters_data->data.push_back(Parameter::Data(Parameter::Template { QualifiedName(name), min, central, max, latex, unit }, idx));
                        parameters_map[name] = idx;
                        parameters.push_back(Parameter(parameters_data, idx));
                    }
                }
            }
            catch (std::exception & e)
            {
                throw ParameterInputFileParseError(file, e.what());
            }
        }
    };

    Parameters::Parameters(Implementation<Parameters> * imp) :
        PrivateImplementationPattern<Parameters>(imp)
    {
    }

    Parameters::~Parameters()
    {
    }

    Parameters
    Parameters::clone() const
    {
        return Parameters(new Implementation<Parameters>(*_imp));
    }

    Parameter
    Parameters::operator[] (const QualifiedName & name) const
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        return Parameter(_imp->parameters_data, i->second);
    }

    Parameter
    Parameters::operator[] (const Parameter::Id & id) const
    {
        if (id >= _imp->parameters.size())
            throw InternalError("Parameters::operator[] (Parameter::Id): invalid id '" + stringify(id) + "'");

        return _imp->parameters[id];
    }

    Parameter::Id
    Parameters::declare(const QualifiedName & name, const std::string & latex, Unit unit, const double & value, const double & min, const double & max)
    {
        return ParameterDefaults::instance()->declare(name, Parameter::Template { name, min, value, max, latex, unit });
    }

    Parameter
    Parameters::declare_and_insert(const QualifiedName & name, const std::string & latex, Unit unit, const double & value, const double & min, const double & max)
    {
        // is the parameter already available?
        auto i(_imp->parameters_map.find(name));
        if (_imp->parameters_map.end() != i)
        {
            Log::instance()->message("[parameters.declare_and_insert]", ll_error)
                << "Parameter '" << name << "' is already declared, returning existing instance; check your code for conflicting duplicate declarations";

            return Parameter(_imp->parameters_data, i->second);
        }

        // declare the new parameter ...
        ParameterDefaults::instance()->declare(name, Parameter::Template { name, min, value, max, latex, unit });

        // ... and insert it into this parameter set ...
        unsigned idx = _imp->parameters.size();
        _imp->parameters_data->data.push_back(Parameter::Data(Parameter::Template { name, min, value, max, latex, unit }, idx));
        _imp->parameters_map[name] = idx;
        _imp->parameters.push_back(Parameter(_imp->parameters_data, idx));

        // ... before returning it
        return _imp->parameters.back();
    }

    void
    Parameters::redirect(const QualifiedName & name, const Parameter::Id & id)
    {
        ParameterDefaults::instance()->redirect(name, id);
    }

    void
    Parameters::redirect_and_apply(const QualifiedName & name, const Parameter::Id & id)
    {
        // redirect the parameter to the new id ...
        ParameterDefaults::instance()->redirect(name, id);

        // ... and apply the change to this parameter set ...
        auto i = _imp->parameters_map.find(name);
        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        i->second = id;
    }

    void
    Parameters::set(const QualifiedName & name, const double & value)
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        _imp->parameters_data->data[i->second].value = value;
    }

    bool
    Parameters::has(const QualifiedName & name)
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            return false;
        else return true;
    }

    Parameters::Iterator
    Parameters::begin() const
    {
        return Parameters::Iterator(_imp->parameters.begin());
    }

    Parameters::Iterator
    Parameters::end() const
    {
        return Parameters::Iterator(_imp->parameters.end());
    }

    Parameters::SectionIterator
    Parameters::begin_sections() const
    {
        return ParameterDefaults::instance()->sections().begin();
    }

    Parameters::SectionIterator
    Parameters::end_sections() const
    {
        return ParameterDefaults::instance()->sections().end();
    }

    bool
    Parameters::operator!= (const Parameters & rhs) const
    {
        return rhs._imp.get() != this->_imp.get();
    }

    Parameters
    Parameters::Defaults()
    {
        auto imp = new Implementation<Parameters>();

        return Parameters(imp);
    }

    void
    Parameters::override_from_file(const std::string & file)
    {
        _imp->override_from_file(file);
    }

    Parameter::Parameter(const std::shared_ptr<Parameters::Data> & parameters_data, unsigned index) :
        _parameters_data(parameters_data),
        _index(index)
    {
    }

    Parameter::Parameter(const Parameter & other) :
        _parameters_data(other._parameters_data),
        _index(other._index)
    {
    }

    Parameter::~Parameter()
    {
    }

    MutablePtr
    Parameter::clone() const
    {
        return MutablePtr(new Parameter(_parameters_data, _index));
    }

    Parameter::operator double () const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::operator() () const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::evaluate() const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::evaluate_generator() const
    {
        return _parameters_data->data[_index].generator_value;
    }

    const Parameter &
    Parameter::operator= (const double & value)
    {
        _parameters_data->data[_index].value = value;

        return *this;
    }

    void
    Parameter::set(const double & value)
    {
        _parameters_data->data[_index].value = value;
    }

    void
    Parameter::set_generator(const double & value)
    {
        _parameters_data->data[_index].generator_value = value;
    }

    const double &
    Parameter::central() const
    {
        return _parameters_data->data[_index].central;
    }

    const double &
    Parameter::max() const
    {
        return _parameters_data->data[_index].max;
    }

    void
    Parameter::set_max(const double & value)
    {
        _parameters_data->data[_index].max = value;
    }

    const double &
    Parameter::min() const
    {
        return _parameters_data->data[_index].min;
    }

    void
    Parameter::set_min(const double & value)
    {
        _parameters_data->data[_index].min = value;
    }

    const std::string &
    Parameter::name() const
    {
        return _parameters_data->data[_index].name.str();
    }

    const std::string &
    Parameter::latex() const
    {
        return _parameters_data->data[_index].latex;
    }

    Unit
    Parameter::unit() const
    {
        return _parameters_data->data[_index].unit;
    }

    Parameter::Id
    Parameter::id() const
    {
        return _parameters_data->data[_index].id;
    }

    /* ParameterUser */

    template <>
    struct WrappedForwardIteratorTraits<ParameterUser::ConstIteratorTag>
    {
        using UnderlyingIterator = std::set<Parameter::Id>::const_iterator;
    };
    template class WrappedForwardIterator<ParameterUser::ConstIteratorTag, const Parameter::Id>;

    ParameterUser::ConstIterator
    ParameterUser::begin() const
    {
        return ConstIterator(_ids.cbegin());
    }

    ParameterUser::ConstIterator
    ParameterUser::end() const
    {
        return ConstIterator(_ids.cend());
    }

    void
    ParameterUser::drop(const Parameter::Id & id)
    {
        _ids.erase(id);
    }

    void
    ParameterUser::uses(const Parameter::Id & id)
    {
        _ids.insert(id);
    }

    void
    ParameterUser::uses(const ParameterUser & other)
    {
        _ids.insert(other._ids.cbegin(), other._ids.cend());
    }

    UsedParameter::UsedParameter(const Parameter & parameter, ParameterUser & user) :
        Parameter(parameter)
    {
        user.uses(parameter.id());
    }

    UnknownParameterError::UnknownParameterError(const QualifiedName & name) throw () :
        Exception("Unknown parameter: '" + name.full() + "'")
    {
    }

    ParameterInputFileParseError::ParameterInputFileParseError(const std::string & file, const std::string & msg) throw () :
        Exception("Malformed parameter input file '" + file + "': " + msg)
    {
    }

    ParameterInputFileNodeError::ParameterInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw () :
        Exception("Malformed parameter input file '" + file + "': Node '" + node + "' " + msg)
    {
    }

    ParameterInputDuplicateError::ParameterInputDuplicateError(const std::string & file, const std::string & node) throw () :
        Exception("Malformed parameter input file '" + file + "': Duplicate entry for parameter '" + node + "'")
    {
    }
}
