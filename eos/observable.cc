/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2017, 2021-2022 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2018, 2019 Nico Gubernari
 * Copyright (c) 2024 Lorenz Gärtner
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

#include <eos/b-decays/observables.hh>
#include <eos/c-decays/observables.hh>
#include <eos/form-factors/observables.hh>
#include <eos/meson-mixing/observables.hh>
#include <eos/nonlocal-form-factors/observables.hh>
#include <eos/observable-impl.hh>
#include <eos/observable.hh>
#include <eos/rare-b-decays/observables.hh>
#include <eos/scattering/observables.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/expression-parser-impl.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/observable_stub.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <map>

namespace eos
{
    namespace impl
    {
        std::map<QualifiedName, ObservableEntryPtr> observable_entries;
    }

    ObservableEntries::ObservableEntries() :
        _entries(&impl::observable_entries)
    {
        std::vector<std::function<ObservableSection()>> section_makers = { make_form_factors_section, make_nonlocal_form_factors_section, make_b_decays_section,
                                                                           make_c_decays_section,     make_rare_b_decays_section,         make_meson_mixing_section,
                                                                           make_scattering_section };

        for (const auto & section_maker : section_makers)
        {
            for (const auto & group : section_maker())
            {
                _entries->insert(group.begin(), group.end());
            }
        }
    }

    ObservableEntries::~ObservableEntries() = default;

    void
    ObservableEntries::insert_or_assign(const QualifiedName & key, const std::shared_ptr<const ObservableEntry> & value)
    {
        auto result = _entries->insert_or_assign(key, value);

        if (! result.second)
        {
            Log::instance()->message("[ObservableEntries.insert_or_assign]", ll_warning) << "Entry for observable " << key.str() << " has been replaced.";
        }
    }

    ObservablePtr
    Observable::make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        const std::map<QualifiedName, std::shared_ptr<const ObservableEntry>> & observable_entries = ObservableEntries::instance()->entries();

        // check if 'name' matches a simple observable
        {
            auto i = observable_entries.find(name);
            if (observable_entries.end() != i)
            {
                return i->second->make(parameters, kinematics, name.options() + _options);
            }
        }

        // check if 'name' matches a parameter
        if (name.options().empty())
        {
            auto i = std::find_if(parameters.begin(), parameters.end(), [&](const Parameter & p) { return p.name() == name.str(); });
            if (parameters.end() != i)
            {
                return ObservablePtr(new ObservableStub(parameters, name));
            }
        }

        throw UnknownObservableError("Expression '" + name.full() + "' is neither a known Observable nor a Parameter");

        return ObservablePtr();
    }

    class ObservableSections : public InstantiationPolicy<ObservableSections, Singleton>
    {
        private:
            std::vector<ObservableSection> _sections;

            ObservableSections()
            {
                // ensure that the observable entries have been generated already
                auto entries = std::distance(ObservableEntries::instance()->entries().begin(), ObservableEntries::instance()->entries().end());
                Log::instance()->message("ObservableSections::ObservableSections()", ll_debug) << "Total number of registered observables: " << entries;

                _sections = std::vector<ObservableSection>({ make_b_decays_section(),
                                                             make_c_decays_section(),
                                                             make_rare_b_decays_section(),
                                                             make_meson_mixing_section(),
                                                             make_nonlocal_form_factors_section(),
                                                             make_form_factors_section(),
                                                             make_scattering_section() });
            }

            ~ObservableSections() = default;

        public:
            friend class InstantiationPolicy<ObservableSections, Singleton>;

            const std::vector<ObservableSection> &
            sections() const
            {
                return _sections;
            }
    };

    /* ObservableEntry */

    ObservableEntry::ObservableEntry() {}

    ObservableEntry::~ObservableEntry() {}

    /* ObservableGroup */

    template <> struct WrappedForwardIteratorTraits<ObservableGroup::ObservableIteratorTag>
    {
            using UnderlyingIterator = std::map<QualifiedName, ObservableEntryPtr>::const_iterator;
    };
    template class WrappedForwardIterator<ObservableGroup::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;

    ObservableGroup::ObservableGroup(Implementation<ObservableGroup> * imp) :
        PrivateImplementationPattern<ObservableGroup>(imp)
    {
    }

    ObservableGroup::~ObservableGroup() = default;

    ObservableGroup::ObservableIterator
    ObservableGroup::begin() const
    {
        return _imp->entries.begin();
    }

    ObservableGroup::ObservableIterator
    ObservableGroup::end() const
    {
        return _imp->entries.end();
    }

    const std::string &
    ObservableGroup::name() const
    {
        return _imp->name;
    }

    const std::string &
    ObservableGroup::description() const
    {
        return _imp->description;
    }

    /* ObservableSection */

    template <> struct WrappedForwardIteratorTraits<ObservableSection::GroupIteratorTag>
    {
            using UnderlyingIterator = std::vector<ObservableGroup>::const_iterator;
    };
    template class WrappedForwardIterator<ObservableSection::GroupIteratorTag, const ObservableGroup &>;

    ObservableSection::ObservableSection(Implementation<ObservableSection> * imp) :
        PrivateImplementationPattern<ObservableSection>(imp)
    {
    }

    ObservableSection::~ObservableSection() = default;

    ObservableSection::GroupIterator
    ObservableSection::begin() const
    {
        return _imp->groups.begin();
    }

    ObservableSection::GroupIterator
    ObservableSection::end() const
    {
        return _imp->groups.end();
    }

    const std::string &
    ObservableSection::name() const
    {
        return _imp->name;
    }

    const std::string &
    ObservableSection::description() const
    {
        return _imp->description;
    }

    /* Observables */

    template <> struct WrappedForwardIteratorTraits<Observables::ObservableIteratorTag>
    {
            using UnderlyingIterator = std::map<QualifiedName, ObservableEntryPtr>::const_iterator;
    };
    template class WrappedForwardIterator<Observables::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;

    template <> struct WrappedForwardIteratorTraits<Observables::SectionIteratorTag>
    {
            using UnderlyingIterator = std::vector<ObservableSection>::const_iterator;
    };
    template class WrappedForwardIterator<Observables::SectionIteratorTag, const ObservableSection &>;

    template <> struct Implementation<Observables>
    {
            std::vector<ObservableSection> observable_sections;

            std::map<QualifiedName, ObservableEntryPtr> observable_entries;

            Implementation() :
                observable_sections(ObservableSections::instance()->sections()),
                observable_entries(ObservableEntries::instance()->entries())
            {
            }
    };

    Observables::Observables() :
        PrivateImplementationPattern<Observables>(new Implementation<Observables>())
    {
    }

    Observables::~Observables() {}

    ObservableEntryPtr
    Observables::operator[] (const QualifiedName & qn) const
    {
        const auto & observable_entries = ObservableEntries::instance()->entries();

        auto i = observable_entries.find(qn);
        if (i != observable_entries.end())
        {
            return i->second;
        }

        throw UnknownObservableError("'" + qn.full() + "' not known");
    }

    Observables::ObservableIterator
    Observables::begin() const
    {
        const auto & observable_entries = ObservableEntries::instance()->entries();

        return ObservableIterator(observable_entries.begin());
    }

    Observables::ObservableIterator
    Observables::end() const
    {
        const auto & observable_entries = ObservableEntries::instance()->entries();

        return ObservableIterator(observable_entries.end());
    }

    Observables::SectionIterator
    Observables::begin_sections() const
    {
        return SectionIterator(_imp->observable_sections.begin());
    }

    Observables::SectionIterator
    Observables::end_sections() const
    {
        return SectionIterator(_imp->observable_sections.end());
    }

    void
    Observables::insert(const QualifiedName & name, const std::string & latex, const Unit & unit, const Options & forced_options, const std::string & input) const
    {
        eos::exp::Expression expression;

        using It = std::string::const_iterator;
        ExpressionParser<It> parser;

        It   first(input.begin()), last(input.end());
        bool completed = qi::phrase_parse(first, last, parser, ascii::space, expression) && (first == last);

        if ((! completed) || expression.empty())
        {
            throw ParsingError("Could not parse expression '" + input + "'");
        }

        ExpressionObservableEntry * expression_observable_entry = new ExpressionObservableEntry(name, latex, unit, expression, forced_options);

        if (! expression_observable_entry)
        {
            throw InternalError("Could not create expression '" + input + "'");
        }

        ObservableEntries::instance()->insert_or_assign(name, std::shared_ptr<const ObservableEntry>(expression_observable_entry));
    }

    bool
    Observables::has(const QualifiedName & name)
    {
        auto i(_imp->observable_entries.find(name));

        if (_imp->observable_entries.end() == i)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    std::pair<QualifiedName, ObservableEntryPtr>
    make_expression_observable(const char * name, const char * latex, const Unit & unit, const char * _expression)
    {
        using namespace exp;

        const QualifiedName qn(name);
        const std::string   input(_expression);
        Expression          expression;

        {
            bool completed;

            using It = std::string::const_iterator;
            ExpressionParser<It> parser;

            It first(input.begin()), last(input.end());
            completed = qi::phrase_parse(first, last, parser, ascii::space, expression) && (first == last);

            if (! completed)
            {
                throw InternalError("Error when parsing expression " + std::string(name) + " in make_expression_observable");
            }
        }

        auto result = std::make_pair(qn, ObservableEntryPtr(new ExpressionObservableEntry(qn, std::string(latex), unit, expression, Options{})));

        impl::observable_entries.insert(result);

        return result;
    }
} // namespace eos
